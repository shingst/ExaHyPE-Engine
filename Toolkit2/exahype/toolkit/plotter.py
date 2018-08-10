import os

import generator

import helper

class PlotterGenerator(generator.Generator):
	def __init__(self,spec,spec_file,verbose):
		generator.Generator.__init__(self,spec,spec_file,verbose)
		
		self._output_directory += "/"+spec["paths"].get("plotter_subdirectory","").strip() # TODO(Dominic): Overwrite with plotter directory for plotters
		if not os.path.exists(self._output_directory):
			print("Create plotter directory '"+self._output_directory+"'")
			os.mkdirs(self.__output_directory)
		else:
			print("Plotter directory '"+self._output_directory+"' already exists")
	
	def generate_all_plotters(self, solver_num, solver):
		context = {}
		context["project_name"]=self._spec["project_name"]
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
				generator.Generator.render_template(self,template_map[user_defined][file_path],context,file_path,"Generated plotter file",False)
