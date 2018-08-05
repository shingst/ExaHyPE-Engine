import jinja2

class SolverGenerator():
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
			kernel     = solver["kernel"]
			print(solver["kernel"])
			
			context = { }
			context["project"]=spec["project_name"];
			
			context["dimension"]=domain["dimension"];
			
			context["solver" ]          =solver["name"];
			context["abstractSolver" ]  ="Abstract"+solver["name"];
			#context["linearOrNonlinear"]= , isLinear? "Linear" : "Nonlinear");
			#context["language"]         = isFortran? "fortran" : "c");
			#context["order"]            = solver["order"]); # solver specific

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
