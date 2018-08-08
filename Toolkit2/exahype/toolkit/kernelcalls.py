import os
import jinja2

import helper

import plotter

class KernelCallsGenerator():
	_exahype_root     = os.path.dirname(os.path.abspath(__file__+"/../../..")) # adjust when file is moved
	
	_spec_file        = None
	_spec             = None
	_output_directory = None
	_dimensions       = None
	_verbose          = False
	_jinja2_env       = None
	
	def __init__(self,spec,spec_file,verbose):
		self._spec_file            = str(spec_file)
		self._spec                 = spec
		self._output_directory     = self._exahype_root+"/"+spec["paths"]["output_directory"]
		self._plotter_subdirectory = spec["paths"].get("plotter_subdirectory","").strip()
		self._dimensions           = spec["computational_domain"]["dimension"]
		self._verbose              = verbose
		self._jinja2_env           = jinja2.Environment(loader=jinja2.FileSystemLoader(self._exahype_root+"/Toolkit2/exahype/toolkit/templates"),trim_blocks=True)
		
	def write_files(self,context):
		template_map = {
			"KerneCalls.cpp" : "KernelCallsImplementation.template"
		}
		for file_path in template_map:
			try:
				template = self._jinja2_env.get_template(template_map[file_path])
				rendered_output = template.render(data=context)
				with open(self._output_directory+"/"+file_path,"w") as file_handle:
					file_handle.write(rendered_output)
				print("Generated solver registration file '"+file_path+"'")
			except Exception as e:
				raise helper.TemplateNotFound(template_map[file_path])
	
	def generate_solver_registration(self):
		"""
		Write the solver registration (KernelCalls.cpp).
		"""
		for solver in self._spec.get("solvers",[]):
			solver["header_path"] = solver["name"]+".h"
			for plotter in solver.get("plotters",[]):
				plotter["header_path"] = ( self._plotter_subdirectory+"/" + plotter["name"]+".h" ).strip("/")
				plotter["type_as_str"] = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
		
		context = {}
		context["spec_file"]=self._spec_file
		context["subPaths"] =[]
		# todo(JM) optimised kernels subPaths
		# todo(JM) profiler 
		# todo(Sven) solver version
		# todo(Sven) serialised spec file compiled into KernelCalls.cp
		context["solvers"]      = self._spec.get("solvers",[])
		context["project_name"] = self._spec["project_name"]
		self.write_files(context)

