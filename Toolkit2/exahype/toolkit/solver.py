import os
import jinja2

import helper

import plotter

class SolverGenerator():
	"""
	Generates the solver files
	"""
	_exahype_root     = os.path.dirname(os.path.abspath(__file__+"/../../..")) # adjust when file is moved
	
	_project_name      = "unknown"
	_output_directory  = "invalid"
	_plotter_directory = "invalid"
	_dimensions        = -1
	_verbose           = False
	_jinja2_env        = None
	
	def __init__(self,spec,verbose):
		self._project_name      = spec["project_name"]
		self._output_directory  = self._exahype_root+"/"+spec["paths"]["output_directory"]
		self._plotter_directory = self._output_directory+"/"+spec["paths"]["plotter_subdirectory"]
		self._dimensions        = spec["computational_domain"]["dimension"]
		self._verbose           = verbose
		self._jinja2_env        = jinja2.Environment(loader=jinja2.FileSystemLoader(self._exahype_root+"/Toolkit2/exahype/toolkit/templates"),trim_blocks=True)
		
		if not os.path.exists(self._output_directory):
			print("Created output directory '"+self._output_directory+"'")
			os.mkdir(self._output_directory)
		else:
			print("Output directory '"+self._output_directory+"' already exists.")
		
		if not os.path.exists(self._plotter_directory):
			print("Create plotter directory '"+self._plotter_directory+"'")
			os.mkdirs(self._plotter_directory)
		else:
			print("Plotter directory '"+self._plotter_directory+"' already exists")
	
	def write_solver_files(self,solver_map,abstract_solver_map,implementation,context):
		for file_path in solver_map.get(implementation,{}):
			if not os.path.exists(self._output_directory+"/"+file_path):
				try:
					template = self._jinja2_env.get_template(solver_map[file_path])
					rendered_output = template.render(data=context)
					with open(self._output_directory+"/"+file_path,"w") as file_handle:
						file_handle.write(rendered_output)
					print("Generated user solver file '"+file_path+"'")
				except Exception as e:
					raise helper.TemplateNotFound(solver_map[file_path])
			else:
				print("User solver file '"+file_path+"' already exists. Is not overwritten.")
		for file_path in abstract_solver_map.get(implementation,{}):
			try:
				template = self._jinja2_env.get_template(abstract_solver_map[implementation][file_path])
				rendered_output = template.render(data=context)
				with open(self._output_directory+"/"+file_path,"w") as file_handle:
					file_handle.write(rendered_output)
				print("Generated abstract solver file '"+file_path+"'")
			except Exception as e:
				raise helper.TemplateNotFound(solver_map[file_path])
	
	def create_solver_context(self,solver):
		context = {}
		context["project"]        = self._project_name
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
	
	def generate_aderdg_solver_files(self,solver,dmp_observables):
		"""
		Generate user solver and abstract solver header and source files for an ADER-DG solver.
		
		Does not overwrite user solver files if they already exist.
		"""
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

	def generate_fv_solver_files(self,solver,patch_size):
		"""
		Generate user solver and abstract solver header and source files for an FV solver.
		
		Does not overwrite user solver files if they already exist.
		"""
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
		print(solver_map)
		
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
	
	def generate_limiting_aderdg_solver_files(self,solver):
		"""
		Generate user solver and abstract solver header and source files for an Limiting-ADER-DG solver.
		Further generate those files for the wrapped ADER-DG and FV solver coupled by 
		the Limiting-ADER-DG solver. 
		
		Does not overwrite user solver files if they already exist.
		"""
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
	
	def generate_all_plotters(self, solver_num, solver):
		context = {}
		context["project_name"]=self._project_name
		context["solver"]=solver
		
		plotters = solver["plotters"]
		for j,plotter in enumerate(solver.get("plotters",[])):
			print("Generating plotter[%d] = %s for solver" % (j, plotter["name"]))
			plotter["number_of_variables"]=helper.count_variables(helper.parse_variables(plotter,"variables"))
			
			context["plotter"]=plotter
			
			type_as_str  = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
			user_defined = type_as_str=="user::defined"
			
			template_map = {
				True : {	plotter["name"]+".h" : "plotters/UserDefinedDeviceHeader.template",
									plotter["name"]+".cpp" : "plotters/UserDefinedDeviceImplementation.template" },
				False: {	plotter["name"]+".h" : "plotters/UserOnTheFlyPostProcessingHeader.template",
									plotter["name"]+".cpp" : "plotters/UserOnTheFlyPostProcessingImplementation.template" }
			}
			
			for file_path in template_map.get(user_defined,{}):
				if not os.path.exists(self._plotter_directory+"/"+file_path):
					try:
						template = self._jinja2_env.get_template(template_map[user_defined][file_path])
						rendered_output = template.render(data=context)
						with open(self._plotter_directory+"/"+file_path,"w") as file_handle:
							file_handle.write(rendered_output)
						print("Generated plotter file '"+file_path+"'")
					except Exception as e:
						raise helper.TemplateNotFound(template_map[user_defined][file_path])
				else:
					print("Plotter file '"+file_path+"' already exists. Is not overwritten.")
	
	def generate_all_solvers(self,spec):
		"""
		Generate user and abstract solver files for all solvers found in the
		specification.
		"""
		for i, solver in enumerate(spec.get("solvers",[])):
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
			
			p = self.generate_all_plotters(i, solver)
			
			
	def generate_solver_registration(self,spec):
		"""
		Write the solver registration (KernelCalls.cpp).
		"""
		for i, solver in enumerate(spec.get("solvers",[])):
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
