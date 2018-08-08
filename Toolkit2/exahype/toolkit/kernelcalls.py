import os
import jinja2

import helper

import plotter

def KernelCallsGenerator(a,b):
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
	
	def generate_solver_registration(self,spec):
		"""
		Write the solver registration (KernelCalls.cpp).
		"""
		for i, solver in enumerate(spec.get("solvers",[])):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
		
			if solver["type"]=="ADER-DG":
				self.generate_aderdg_solver_files(
					solver=solver,\
					solver_name=solver["name"],\
					dmp_observables=0)
			elif solver["type"]=="Finite-Volumes":
				self.generate_fv_solver_files(
					solver=solver,\
					solver_name=solver["name"],\
					patch_size=solver["patch_size"])
			elif solver["type"]=="Limiting-ADER-DG":
				self.generate_limiting_aderdg_solver_files(solver)
