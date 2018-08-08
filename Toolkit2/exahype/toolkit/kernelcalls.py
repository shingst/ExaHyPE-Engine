import os
import jinja2

import helper

import plotter

class KernelCallsGenerator():
	_exahype_root     = os.path.dirname(os.path.abspath(__file__+"/../../..")) # adjust when file is moved
	
	_project_name     = "unknown"
	_output_directory = "invalid"
	_dimensions       = -1
	_verbose          = False
	_jinja2_env       = None
	
	def __init__(self,spec,verbose):
		self._project_name      = spec["project_name"]
		self._output_directory  = self._exahype_root+"/"+spec["paths"]["output_directory"]
		self._plotter_directory = self._output_directory+spec["paths"].get("plotter_subdirectory","")
		self._dimensions        = spec["computational_domain"]["dimension"]
		self._verbose           = verbose
		self._jinja2_env        = jinja2.Environment(loader=jinja2.FileSystemLoader(self._exahype_root+"/Toolkit2/exahype/toolkit/templates"),trim_blocks=True)
	
	def write_files(self,context):
		template_map = {
			"KerneCalls.cpp" : "KernelCallsImplementation.template"
		}
		for file_path in template_map:
			try:
				template = self._jinja2_env.get_template(template_map[file_path])
				rendered_output = template.render(context)
				with open(self._output_directory+"/"+file_path,"w") as file_handle:
					file_handle.write(rendered_output)
				print("Generated solver registration file '"+file_path+"'")
			except Exception as e:
				raise helper.TemplateNotFound(template_map[file_path])
	
	def generate_solver_registration(self,spec):
		"""
		Write the solver registration (KernelCalls.cpp).
		"""
		context = {}
		context["subPaths"]=0
		# todo optimised kernels subPaths
		
		for solver in spec.get("solvers",[]):
			solver["header_path"] = self._output_directory+"/"+solver["name"]+".h"
			for plotter in solver.get("plotters",[]):
				plotter["header_path"] = self._plotter_directory+"/"+plotter["name"]+".h"
		
		print(spec)
		
	def write_files(self,context):
		pass

