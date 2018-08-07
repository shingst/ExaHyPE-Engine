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

def convert_to_dict(list_of_str_or_dict):
	"""
	Example:
	
	Converts '[ 'generic', 'usestack', { 'maxpicarditer' : 1 } ]' to
	'{ 'generic' : True , 'usestack' : True , 'maxpicarditer' : 1 }'
	"""
	result = {}
	for i,item in enumerate(list_of_str_or_dict):
		if type(item) is str:
			result[item] = True
		elif type(item) is dict:
			result[list(item.keys())[0]] = list(item.values())[0]
	return result

class SolverGenerator():
	
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
	
	def create_kernel_context(self,kernel_optimisations,kernel_terms,template_bool_map,template_int_map):
		"""
		Helper function.
		"""
		context = {}
		for key in template_bool_map:
			context[template_bool_map[key]] = \
				key in kernel_optimisations or key in kernel_terms
		for key in template_int_map:
			context[template_int_map[key]] = \
				max(kernel_optimisations.get(key,0),kernel_terms.get(key,0))
		return context
		
	def create_aderdg_kernel_context(self,kernel):
		kernel_terms         = convert_to_dict(kernel["terms"])
		kernel_optimisations = convert_to_dict(kernel["optimisations"])
		
		template_bool_map = {
			"linear"             : "isLinear",
			"nonlinear"          : "isNonlinear",
			"Fortran"            : "isFortran",
			"flux"               : "useFlux",
			"source"             : "useSource",
			"ncp"                : "useNCP",
			"pointsources"       : "usePointSources",
			"materialparameters" : "useMaterialParam",
			"notimeavg"          : "noTimeAveraging",
			"patchwiseadjust"    : "patchwiseAdjust",
			"usestack"           : "tempVarsOnStack",
			"maxpicarditer"      : "useMaxPicardIterations",
			"flops"              : "countFlops"
		}
		
		template_int_map = {
			"pointsources"    : "numberOfPointSources",
			"maxpicarditer"   : "maxPicardIterations"
		}
		
		context = self.create_kernel_context(kernel_optimisations,kernel_terms,template_bool_map,template_int_map)
		
		print(" ".join([str(x) for x in range(0,10)]))
		
		context[template_bool_map[kernel["type"]]] = True
		context["isFortran"]                       = True if kernel["language"] is "Fortran" else False;
		return context
		
	def create_fv_kernel_context(self,kernel):
		kernel_type          = kernel["type"]
		kernel_terms         = convert_to_dict(kernel["terms"])
		kernel_optimisations = convert_to_dict(kernel["optimisations"])
		
		template_bool_map = {
			"flux"               : "useFlux",
			"source"             : "useSource",
			"ncp"                : "useNCP",
			"pointsources"       : "usePointSources",
			"materialparameters" : "useMaterialParam", # todo not implemented yet
			"patchwiseadjust"    : "patchwiseAdjust",  # todo not implemented yet
			"usestack"           : "tempVarsOnStack",
			"flops"              : "countFlops"        # todo not implemented yet
		}
		
		template_int_map = {
			"pointsources"    : "numberOfPointSources"
		}
		
		context = self.create_kernel_context(kernel_optimisations,kernel_terms,template_bool_map,template_int_map)
		return context
	
	def create_solver_context(self,solver,solverName):
		context = {}
		context["project"]        =self._project_name
		context["dimensions"]     =self._dimensions
		context["solver"]         = solverName
		context["abstractSolver"] = "Abstract"+solverName
		
		context["numberOfVariables"]          = numberOfVariables(parseVariables(solver,"variables"));
		context["numberOfMaterialParameters"] = numberOfVariables(parseVariables(solver,"material_parameters"));
		context["numberOfGlobalObservables"]  = numberOfVariables(parseVariables(solver,"global_observables"));
		return context
	
	def generate_aderdg_solver_files(self,solver):
		#aderdg_types=[\
		#	"Legendre",\
		#	"Lobatto",\
		#	"user"\
		#]
		aderdg_context = self.create_solver_context(solver,solver["name"])
		aderdg_context.update(self.create_aderdg_kernel_context(solver["aderdg_kernel"]))
		print(aderdg_context)

	def generate_fv_solver_files(self,solver):
		fv_context = self.create_solver_context(solver,solver["name"])
		fv_context.update(self.create_fv_kernel_context(solver["fv_kernel"],verbose))
		
	def generate_limiting_aderdg_solver_files(self,solver):
		# aderdg
		aderdg_context = self.create_solver_context(solver,solver["name"]+"_ADERDG")
		aderdg_context.update(self.create_aderdg_kernel_context(solver["aderdg_kernel"]))
		# fv
		fv_context = self.create_solver_context(solver,solver["name"]+"_FV")
		fv_context.update(self.create_fv_kernel_context(solver["fv_kernel"]))
		# limiter
		limiting_context = self.create_solver_context(solver,solver["name"])
	
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
