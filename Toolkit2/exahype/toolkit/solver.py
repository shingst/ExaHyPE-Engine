import os
import jinja2

import helper

import generator
import plotter

##
# Generates the solver files
class SolverGenerator(generator.Generator):
	_plotter_generator = None
	
	def __init__(self,spec,spec_file,verbose):
		generator.Generator.__init__(self,spec,spec_file,verbose)
		self._plotter_generator = plotter.PlotterGenerator(spec,spec_file,verbose)
	
	def create_solver_context(self,solver):
		context = {}
		context["project"]        = self._spec["project_name"]
		context["solver"]         = solver
		
		context["dimensions"] = self._dimensions
		
		nVar          = helper.count_variables(helper.parse_variables(solver,"variables"))
		nParam        = helper.count_variables(helper.parse_variables(solver,"material_parameters"))
		nGlobalObs    = helper.count_variables(helper.parse_variables(solver,"global_observables"))
		nPointSources = len(solver["point_sources"])
		
		context["numberOfVariables"]          = nVar
		context["numberOfMaterialParameters"] = nParam
		context["numberOfGlobalObservables"]  = nGlobalObs
		context["numberOfPointSources"]       = nPointSources
		
		context["range_0_nDim"]          = range(0,self._dimensions)
		context["range_0_nVar"]          = range(0,nVar)
		context["range_0_nVarParam"]     = range(0,nVar+nParam)
		context["range_0_nGlobalObs"]    = range(0,nGlobalObs)    # nGlobalObs might be 0
		context["range_0_nPointSources"] = range(0,nPointSources) # nPointSources might be 0
		
		return context
	
	def write_solver_files(self,solver_map,abstract_solver_map,implementation,context):
		for file_path in solver_map.get(implementation,{}):
			generator.Generator.render_template(self,solver_map[implementation][file_path],context,file_path,"Generated user solver file",False)
		for file_path in abstract_solver_map.get(implementation,{}):
			generator.Generator.render_template(self,abstract_solver_map[implementation][file_path],context,file_path,"Generated abstract solver file",True)
	
	##
	#Generate user solver and abstract solver header and source files for an ADER-DG solver.
	#
	# @note Does not overwrite user solver files if they already exist.
	def generate_aderdg_solver_files(self,solver,dmp_observables):
		solver_name          = solver["name"]+ ( "_ADERDG" if solver["type"]=="Limiting-ADER-DG" else "" )
		abstract_solver_name = "Abstract"+solver_name
		aderdg_context = self.create_solver_context(solver)
		aderdg_context["numberOfDMPObservables"]=dmp_observables;
		
		solver_map = {
			"user"    :  { solver_name+".h"               : "solvers/MinimalADERDGSolverHeader.template", 
			               solver_name+".cpp"             : "solvers/EmptyADERDGSolverImplementation.template" },
			"generic" : {  solver_name+".h"               : "solvers/ADERDGSolverHeader.template", 
			               solver_name+".cpp"             : "solvers/ADERDGSolverInCUserCode.template"},
		}
		solver_map["optimised"] = solver_map["generic"]
		
		abstract_solver_map  = { 
			"user"      :  { },
			"generic"   :  { abstract_solver_name+".cpp"    : "solvers/AbstractGenericADERDGSolverImplementation.template", 
			                 abstract_solver_name+".h"      : "solvers/AbstractGenericADERDGSolverHeader.template" },
			"optimised" :  { abstract_solver_name+".cpp"    : "solvers/AbstractOptimisedADERDGSolverImplementation.template", 
			                 abstract_solver_name+".h"      : "solvers/AbstractOptimisedADERDGSolverHeader.template" }
		}
		
		implementation = solver["aderdg_kernel"].get("implementation","generic")
		try:
			self.write_solver_files(solver_map,abstract_solver_map,implementation,aderdg_context)
		except Exception:
			raise
		
		if implementation=="optimised":
			print("ERROR: not implemented yet '"+file_path+"'",file=sys.stderr)
			sys.exit(-1)

	##
	# Generate user solver and abstract solver header and source files for an FV solver.
	#
	# @note Does not overwrite user solver files if they already exist.
	#
	def generate_fv_solver_files(self,solver,patch_size):
		solver_name          = solver["name"]+ ( "_FV" if solver["type"]=="Limiting-ADER-DG" else "" )
		abstract_solver_name = "Abstract"+solver_name
		fv_context = self.create_solver_context(solver)
		ghost_layer_width = { "godunov" : 1, "musclhancock" : 2 }
		fv_context["ghostLayerWidth"]=ghost_layer_width[solver["fv_kernel"]["scheme"]]
		
		solver_map = {
			"user"    : { solver_name+".h"   : "solvers/MinimalFiniteVolumesSolverHeader.template", 
			              solver_name+".cpp" : "solvers/EmptyFiniteVolumesSolverImplementation.template" },
			"generic" : { solver_name+".h"   : "solvers/FiniteVolumesHeader.template", 
			              solver_name+".cpp" : "solvers/FiniteVolumesInCUserCode.template"},
		}
		
		abstract_solver_map  = { 
			"generic"   :  { abstract_solver_name+".cpp" : "solvers/AbstractGenericFiniteVolumesSolverImplementation.template", 
			                 abstract_solver_name+".h"   : "solvers/AbstractGenericFiniteVolumesSolverHeader.template" }
		}
		
		implementation = solver["fv_kernel"].get("implementation","generic")
		try:
			self.write_solver_files(solver_map,abstract_solver_map,implementation,fv_context)
		except Exception:
			raise
		
		if implementation=="optimised":
			print("ERROR: not implemented yet '"+file_path+"'",file=sys.stderr)
			sys.exit(-1)
	##
	# Generate user solver and abstract solver header and source files for an Limiting-ADER-DG solver.
	# Further generate those files for the wrapped ADER-DG and FV solver coupled by 
	# the Limiting-ADER-DG solver. 
	#
	# @note Does not overwrite user solver files if they already exist.
	def generate_limiting_aderdg_solver_files(self,solver):
		solver_name          = solver["name"]
		abstract_solver_name = "Abstract"+solver_name
		# aderdg
		self.generate_aderdg_solver_files(
			solver=solver,\
			dmp_observables=solver["limiter"].get("dmp_observables",0))
		# fv
		self.generate_fv_solver_files(
			solver=solver,\
			patch_size=2*solver["order"]+1)
		# limiter
		limiter_context = self.create_solver_context(solver)
		
		solver_map = {  }
		abstract_solver_map  = { 
			"generic"   :  { abstract_solver_name+".cpp" : "solvers/AbstractGenericLimiterSolverImplementation.template", 
			                 abstract_solver_name+".h"   : "solvers/AbstractGenericLimiterSolverHeader.template" },
			"optimised" :  { abstract_solver_name+".cpp" : "solvers/AbstractOptimisedLimiterSolverImplementation.template", 
			                 abstract_solver_name+".h"   : "solvers/AbstractOptimisedLimiterSolverHeader.template" }
		}
		implementation = solver["fv_kernel"].get("implementation","generic")
		try:
			self.write_solver_files(solver_map,abstract_solver_map,implementation,limiter_context)
		except Exception:
			raise
	##
	# Generate user and abstract solver files for all solvers found in the
	# specification.
	#
	def generate_all_solvers(self):
		for i, solver in enumerate(self._spec.get("solvers",[])):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
			
			if solver["type"]=="ADER-DG":
				self.generate_aderdg_solver_files(
					solver=solver,\
					dmp_observables=0)
			elif solver["type"]=="Finite-Volumes":
				self.generate_fv_solver_files(
					solver=solver,\
					patch_size=solver["patch_size"])
			elif solver["type"]=="Limiting-ADER-DG":
				self.generate_limiting_aderdg_solver_files(solver)
			
			self._plotter_generator.generate_all_plotters(i, solver)
