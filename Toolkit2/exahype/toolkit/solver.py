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
	
	def __init__(self,spec,verbose):
		project_name = spec["project_name"]
		solvers = spec["solvers"]
		domain  = spec["computational_domain"]
		
		for i, solver in enumerate(solvers):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
			
			solverType = solver["type"]
			
			kernels             = solver["aderdg_kernel"]
			kernel_type         = kernels["type"]
			kernel_terms        = kernels["terms"]
			kernel_optimisation = kernels["optimisations"]
			print(solver["aderdg_kernel"])
			
			context = { }
			context["project"]=spec["project_name"]
			
			context["dimension"]=domain["dimension"]
			
			context["solver" ]           = solver["name"]
			context["abstractSolver" ]   = "Abstract"+solver["name"]
			context["linearOrNonlinear"] = "Linear"  if ( kernels["type"]=="linear" ) else "Nonlinear"
			context["language"]          = "fortran" if ( kernels["language"]=="Fortran" ) else "c"
			
			context["numberOfVariables"]          = numberOfVariables(parseVariables(solver,"variables"));
			context["numberOfMaterialParameters"] = numberOfVariables(parseVariables(solver,"material_parameters"));
			context["numberOfGlobalObservables"]  = numberOfVariables(parseVariables(solver,"global_observables"));
#			context["numberOfDMPObservables"]     = solver.get("dmp_observables",0);
#			context["numberOfPointSources"]       = solver.get("point_sources",0);

#			context["order"]     = solver.get("order",-1);
#			context["patchSize"] = solver.get("patch_size",-1);
			
			print(context)

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
