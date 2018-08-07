import os
import jinja2

import helper

class SolverGenerator():
	"""
	Generates the solver files
	"""
	_exahype_root     = os.path.dirname(os.path.abspath(__file__+"/../../..")) # adjust when file is moved
	
	_project_name     = "unknown"
	_output_directory = "invalid"
	_dimensions       = -1
	_verbose          = False
	_jinja2_env       = None
	
	def __init__(self,spec,verbose):
		self._project_name     = spec["project_name"]
		self._output_directory = self._exahype_root+"/"+spec["paths"]["output_directory"]
		self._dimensions       = spec["computational_domain"]["dimension"]
		self._verbose          = verbose
		self._jinja2_env       = jinja2.Environment(loader=jinja2.FileSystemLoader(self._exahype_root+"/Toolkit2/exahype/toolkit/templates"),trim_blocks=True)
		
		if not os.path.exists(self._output_directory):
			raise BadSpecificationFile()
	
	def parseVariables(self,solver,field):
		"""
		Parse the 'variables' parameter and similarly
		structured parameters ('material_parameters','global_observables')
		"""
		if field in solver:
			if type(solver[field]) is int:
				return [{"name" : "Q", "multiplicity" : solver[field]}]
			else:
				return solver[field]
		else:
			return [{"name" : "Q", "multiplicity" : 0}]
	
	def numberOfVariables(self,variables):
		"""
		Sum up the multiplicities.
		"""
		number = 0;
		for variable in variables:
			number += variable["multiplicity"]
		return number
	
	def create_kernel_terms_context(self,kernel_terms):
		"""
		Helper function.
		"""
		term_map = {
			"flux"               : "useFlux",
			"source"             : "useSource",
			"ncp"                : "useNCP",
			"pointsources"       : "usePointSources",
			"materialparameters" : "useMaterialParam",
		}
		context = {}
		for key in term_map:
			context[term_map[key]] = key in kernel_terms
		
		return context
	
	def create_aderdg_kernel_context(self,kernel):
		if kernel["type"] != "user":
			context = self.create_kernel_terms_context(kernel["terms"])
			
			context["tempVarsOnStack"]        = kernel.get("allocate_temporary_arrays","heap") is "stack"
			context["patchwiseAdjust"]        = kernel.get("adjust_solution","pointwise")      is "patchwise"
			context["useMaxPicardIterations"] = kernel.get("space_time_predictor",{}).get("maxpicarditer",0)!=0
			context["maxPicardIterations"]    = kernel.get("space_time_predictor",{}).get("maxpicarditer",0)
			context["noTimeAveraging"]        = kernel.get("space_time_predictor",{}).get("notimeavg",False)
			# context["useCERK"]              = kernel.get("space_time_predictor",{}).get("cerkguess",False)
			# context["userConverter"]        = "converter" in kernel.get("optimised_kernel_debugging",[])
			context["countFlops"]             = "flops" in kernel.get("optimised_kernel_debugging",[])
			#context["ADERDGBasis"]           = kernel["basis"]
			
			context["isLinear"]    = kernel.get("type","nonlinear")=="linear"
			context["isNonlinear"] = kernel.get("type","nonlinear")=="nonlinear"
			context["isFortran"]   = True if kernel["language"] is "Fortran" else False;
			return context
		else:
			return {}
		
	def create_fv_kernel_context(self,kernel):
		if kernel["type"] != "user":
			context = self.create_kernel_terms_context(kernel_terms)
			
			context["tempVarsOnStack"]   = kernel.get("allocate_temporary_arrays","heap") is "stack"
			context["patchwiseAdjust"]   = kernel.get("adjust_solution","pointwise")      is "patchwise"
			context["countFlops"]        = "flops" in kernel.get("optimised_kernel_debugging",[])
			
			context["finiteVolumesType"] = kernel["type"]
			return context
		else:
			return {}
	
	def create_solver_context(self,solver,solverName):
		context = {}
		context["project"]        = self._project_name
		context["solver"]         = solverName
		context["abstractSolver"] = "Abstract"+solverName
		
		context["dimensions"] = self._dimensions
		
		nVar          = self.numberOfVariables(self.parseVariables(solver,"variables"))
		nParam        = self.numberOfVariables(self.parseVariables(solver,"material_parameters"))
		nGlobalObs    = self.numberOfVariables(self.parseVariables(solver,"global_observables"))
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
	
	def generate_aderdg_solver_files(self,solver,solver_name,dmp_observables):
		#aderdg_basis = solver["aderdg_kernel"]["basis"]
		abstract_solver_name = "Abstract"+solver_name
		aderdg_context = self.create_solver_context(solver,solver_name)
		aderdg_context.update(self.create_aderdg_kernel_context(solver["aderdg_kernel"]))
		aderdg_context["order"]=solver["order"];
		aderdg_context["numberOfDMPObservables"]=dmp_observables;
		
		aderdg_solver_map = {
			self._output_directory+"/"+solver_name+".h"               : "solvers/ADERDGSolverHeader.template", 
			self._output_directory+"/"+solver_name+".cpp"             : "solvers/ADERDGSolverInCUserCode.template"
		}
		aderdg_abstract_solver_map = { 
			"generic"   :  { self._output_directory+"/"+abstract_solver_name+".cpp"    : "solvers/AbstractGenericADERDGSolverImplementation.template", 
			                 self._output_directory+"/"+abstract_solver_name+".h"      : "solvers/AbstractGenericADERDGSolverHeader.template" },
			"optimised" :  { self._output_directory+"/"+abstract_solver_name+".cpp"    : "solvers/AbstractOptimisedADERDGSolverImplementation.template", 
			                 self._output_directory+"/"+abstract_solver_name+".h"      : "solvers/AbstractOptimisedADERDGSolverHeader.template" }
		}
		
		code = solver["aderdg_kernel"].get("code","generic")
		for file_path in aderdg_solver_map:
			if not os.path.exists(file_path):
				try:
					template = self._jinja2_env.get_template(aderdg_solver_map[file_path])
					rendered_output = template.render(aderdg_context)
					with open(file_path,"w") as file_handle:
						file_handle.write(rendered_output)
					print("Generated user solver file '"+file_path+"'")
				except Exception as e:
					raise helper.TemplateNotFound(aderdg_solver_map[file_path])
			else:
				print("File '"+file_path+"' already exists. Is not overwritten.")
		for file_path in aderdg_abstract_solver_map[code]:
			try:
				template = self._jinja2_env.get_template(aderdg_abstract_solver_map[code][file_path])
				rendered_output = template.render(aderdg_context)
				with open(file_path,"w") as file_handle:
					file_handle.write(rendered_output)
				print("Generated abstract solver file '"+file_path+"'")
			except Exception as e:
				raise helper.TemplateNotFound(aderdg_solver_map[file_path])
		
		if code=="optimised":
			print("ERROR: not implemented yet '"+file_path+"'",file=sys.stderr)
			sys.exit(-1)

	def generate_fv_solver_files(self,solver,solver_name,patch_size):
		abstract_solver_name = "Abstract"+solver_name
		fv_context = self.create_solver_context(solver,solver_name)
		fv_context.update(self.create_fv_kernel_context(solver["fv_kernel"],verbose))
		fv_context["patch_size"]=patch_size
		
	def generate_limiting_aderdg_solver_files(self,solver):
		# aderdg
		self.generate_aderdg_solver_files(
			solver=solver,\
			solver_name=solver["name"]+"_ADERDG",\
			dmp_observables=solver["limiter"].get("dmp_observables",0))
		# fv
		self.generate_fv_solver_files(
			solver=solver,\
			solver_name=solver["name"]+"_FV",\
			patch_size=2*solver["order"]+1)
		# limiter
		limiter_context = self.create_solver_context(solver,solver["name"])
		limiter_context["ADERDGSolver"]         = solver["name"]+"_ADERDG"
		limiter_context["FVSolver"]             = solver["name"]+"_FV"
		limiter_context["ADERDGAbstractSolver"] = "Abstract"+limiter_context["ADERDGSolver"]
		limiter_context["FVAbstractSolver"]     = "Abstract"+limiter_context["FVSolver"]
		
	def generate_all_plotters(self, solver_num, solver):
		plotters = solver["plotters"]
		for j,plotter in enumerate(solver.get("plotters",[])):
			print("Generating plotter[%d] = %s for solver" % (j, plotter["name"]))
		
	def generate_all_solvers(self,spec):
		for i, solver in enumerate(spec.get("solvers",[])):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
			
			if solver["type"]=="ADER-DG":
				self.generate_aderdg_solver_files(
					solver=solver,\
					solver_name=solver["name"],\
					dmp_observables=0)
			elif solver["type"]=="Finite-Volumes":
				self.generate_aderdg_solver_files(
					solver=solver,\
					solver_name=solver["name"],\
					patch_size=solver["patch_size"])
			elif solver["type"]=="Limiting-ADER-DG":
				self.generate_limiting_aderdg_solver_files(solver)
			
			self.generate_all_plotters(i, solver)
