import os
import jinja2

import helper

import generator

class KernelCallsGenerator(generator.Generator):
	def __init__(self,spec,spec_file,verbose):
		generator.Generator.__init__(self,spec,spec_file,verbose)
	
	##
	# Write the solver registration (KernelCalls.cpp).
	def generate_solver_registration(self):
		plotter_subdirectory = self._spec["paths"].get("plotter_subdirectory","").strip()
		
		for solver in self._spec.get("solvers",[]):
			solver["header_path"] = solver["name"]+".h"
			solver["variables_as_str"]            = helper.variables_to_str(solver,"variables")
			solver["material_parameters_as_str"]  = helper.variables_to_str(solver,"material_parameters")
			solver["global_observables_as_str"]   = helper.variables_to_str(solver,"global_observables")
			for plotter in solver.get("plotters",[]):
				plotter["header_path"]      = ( plotter_subdirectory+"/" + plotter["name"]+".h" ).strip("/")
				plotter["type_as_str"]      = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
				plotter["variables_as_str"] = helper.variables_to_str(plotter,"variables")
		
		context = {}
		context["spec_file"]        = self._spec_file
		context["spec_file_as_hex"] = "0x2F" # todo(Sven) do the conversion
		context["subPaths"]         = []
		# todo(JM) optimised kernels subPaths
		# todo(JM) profiler 
		# todo(Sven) serialised spec file compiled into KernelCalls.cp
		context["solvers"]      = self._spec.get("solvers",[])
		context["project_name"] = self._spec["project_name"]
		
		generator.Generator.render_template(self,"KernelCallsImplementation.template",context,"KernelCalls.cpp","Generated solver registration",True)
