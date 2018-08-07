import jinja2

def parseVariables(solver,field):
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

def numberOfVariables(variables):
	"""
	Sum up the multiplicities.
	"""
	number = 0;
	for variable in variables:
		number += variable["multiplicity"]
	return number
	
class SolverGenerator():
	aderdg_template_bool_map = {
		"linear"          : "isLinear",
		"nonlinear"       : "isNonlinear",
		"Fortran"         : "isFortran",
		"flux"            : "useFlux",
		"source"          : "useSource",
		"ncp"             : "useNCP",
		"pointsources"    : "usePointSources",
		"notimeavg"       : "noTimeAveraging",
		"patchwiseadjust" : "patchwiseAdjust",
		"usestack"        : "tempVarsOnStack",
		"maxpicarditer"   : "useMaxPicardIterations",
		"flops"           : "countFlops"
	}
	
	aderdg_types=[\
		"linear",\
		"nonlinear",\
		"Legendre",\
		"Lobatto",\
		"user"\
	]
	
	aderdg_terms=[\
		"flux",\
		"source",\
		"ncp",\
		"pointsources",\
		"materialparameters"\
	]
	
	aderdg_optimisations=[\
		"generic",\
		"optimised",\
		"notimeavg",\
		"patchwiseadjust",\
		"usestack",\
		"maxpicarditer",\
		"fusedsource",\
		"fluxvect",\
		"fusedsourcevect",\
		"cerkguess",\
		"converter",\
		"flops"\
	]
	
	FV_TYPES=[\
		"godunov",\
		"musclhancock",\
		"user",\
	]
	FV_TERMS=[\
		"flux",\
		"source",\
		"ncp",\
		"pointsources",\
		"materialparameters"\
	]
	FV_OPTIMISATIONS=[\
		"generic",\
		"optimised",\
		"patchwiseadjust",\
		"usestack",\
	]
	
	def generate_plotter(self, solver_num, solver):
		plotters = solver["plotters"]
		for j, plotter in enumerate(plotters):
			print("Generating plotter[%d] = %s for solver" % (j, plotter["name"]))
		
			# stub, do something
	
	def create_aderdg_kernel_context(self,kernels):
		kernel_type         = kernels["type"]
		kernel_terms        = kernels["terms"]
		kernel_optimisation = kernels["optimisations"]
		return {}
		
	def create_fv_kernel_context(self,kernels):
		kernel_type         = kernels["type"]
		kernel_terms        = kernels["terms"]
		kernel_optimisation = kernels["optimisations"]
		return {}
	
	def create_solver_context(self,solverName):
		context = {}
		context["project"]        =self._project_name
		context["dimensions"]     =self._dimensions
		context["solver"]         = solverName
		context["abstractSolver"] = "Abstract"+solverName
		return context
	
	def generate_aderdg_solver_files(self,solver):
		aderdg_context = self.create_solver_context(solver["name"])
		aderdg_context.update(self.create_aderdg_kernel_context(solver["aderdg_kernel"]))

	def generate_fv_solver_files(self,solver):
		fv_context = self.create_solver_context(solver["name"])
		fv_context.update(self.create_fv_kernel_context(solver["fv_kernel"],verbose))
		
	def generate_limiting_aderdg_solver_files(self,solver):
		# aderdg
		aderdg_context = self.create_solver_context(solver["name"]+"_ADERDG")
		aderdg_context.update(self.create_aderdg_kernel_context(solver["aderdg_kernel"]))
		# fv
		fv_context = self.create_solver_context(solver["name"]+"_FV")
		fv_context.update(self.create_fv_kernel_context(solver["fv_kernel"]))
		# limiter
		limiting_context = self.create_solver_context(solver["name"])
	
	_project_name = "unknown"
	_dimensions   = -1
	_verbose      = False
	
	def __init__(self,spec,verbose):
		self._project_name = spec["project_name"]
		self._dimensions   = spec["computational_domain"]["dimension"]
		self._verbose      = verbose
	
	def generate_all_solvers(self,spec):
		for i, solver in enumerate(spec["solvers"]):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
			
			if solver["type"]=="ADER-DG":
				self.generate_aderdg_solver_files(solver)
			elif solver["type"]=="Finite-Volumes":
				self.generate_FV_fv_files(solver)
			elif solver["type"]=="Limiting-ADER-DG":
				self.generate_limiting_aderdg_solver_files(solver)
			
#			kernels             = solver["aderdg_kernel"]
#			kernel_type         = kernels["type"]
#			kernel_terms        = kernels["terms"]
#			kernel_optimisation = kernels["optimisations"]
#			print(solver["aderdg_kernel"])
#			
#			
#			context["solver" ]           = solver["name"]
#			context["abstractSolver" ]   = "Abstract"+solver["name"]
#			context["linearOrNonlinear"] = "Linear"  if ( kernels["type"]=="linear" ) else "Nonlinear"
#			context["language"]          = "fortran" if ( kernels["language"]=="Fortran" ) else "c"
#			
#			context["numberOfVariables"]          = numberOfVariables(parseVariables(solver,"variables"));
#			context["numberOfMaterialParameters"] = numberOfVariables(parseVariables(solver,"material_parameters"));
#			context["numberOfGlobalObservables"]  = numberOfVariables(parseVariables(solver,"global_observables"));
#			context["numberOfDMPObservables"]     = solver.get("dmp_observables",0);
#			context["numberOfPointSources"]       = solver.get("point_sources",0);

#			context["order"]     = solver.get("order",-1);
#			context["patchSize"] = solver.get("patch_size",-1);

#//int

#context["numberOfVariables"   , numberOfVariables);
#context["numberOfParameters"  , numberOfParameters);
#context["numberOfObservables" , numberOfObservables);
#context["numberOfPointSources", numberOfPointSources);
#context["maxPicardIterations" , maxPicardIterations);
#
#//boolean
#context["enableProfiler"##, enableProfiler);
#// context["hasConstants"##  , hasConstants);
#context["isLinear"###  , isLinear);
#context["isFortran"### , isFortran);
#context["useFlux"###   , useFlux);
#context["useSource"### , useSource);
#context["useNCP"####, useNCP);
#context["usePointSources"#   , usePointSources);
#context["useMaterialParam"#  , useMaterialParam);
#context["noTimeAveraging"#   , noTimeAveraging);
#context["patchwiseAdjust"#   , patchwiseAdjust);
#context["tempVarsOnStack"#   , tempVarsOnStack);
#context["useMaxPicardIterations", useMaxPicardIterations);
#
#//boolean as String
#context["useFlux_s"## , boolToTemplate(useFlux));
#context["useSource_s"#   , boolToTemplate(useSource));
#context["useNCP_s"##  , boolToTemplate(useNCP));
			
			
			
			# stub, do something
			
			self.generate_plotter(i, solver)
