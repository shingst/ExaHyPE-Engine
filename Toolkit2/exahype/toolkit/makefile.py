import os
import jinja2

import helper

import generator

class MakefileGenerator(generator.Generator):
	def __init__(self,spec,spec_file,verbose):
		generator.Generator.__init__(self,spec,spec_file,verbose)
	
	##
	# Write the application-specific Makefile
	#
	# TODO considered:
	#{{exahypePath}} x
	#{{outputPath}} x
	#{{architecture}}
	#{{alignment}}
	#{{exahypePath}}
	#useOptKernel x
	#useSharedMem x
	#useDistributedMem x
	#useFortran x
	#
	# TODO not considered:
	#useIpcm
	#useLikwid
	#likwidInc
	#
	# todo(JM) optimised kernels subPaths
	# todo(JM) profiler 
	def generate_makefile(self):
		context = {}
		context["project_name"] = self._spec["project_name"]
		context["dimensions"]   = self._dimensions
				
		# dependencies
		context["paths"]        = self._spec.get("paths",[])
		
		context["useSharedMem"]      = "shared_memory" in self._spec;
		context["useDistributedMem"] = "distributed_memory" in self._spec;
		if not "TBB_INC" in os.environ:
			print("WARNING: environment variable TBB_INC not set but required if code is built with TBB");
		if not "TBB_SHLIB" in os.environ:
			print("WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB");
		
		context["useIpcm"]   = False # TODO
		context["useLikwid"] = False # TODO
		context["likwidInc"] = ""    # TODO
		
		# architecture
		context["architecture"] = self._spec.get("architecture","noarch")
		alignment = {
			"snb" : 32, "hsw" : 32, "knc" : 64, "knl" : 64
		}
		context["alignment"] = alignment.get(context["architecture"],16)
		
		# kernels
		useOptKernel = False
		useFortran   = False
		for solver in self._spec["solvers"]:
			if "aderdg_kernel" in solver:
				useOptKernel = useOptKernel or solver["aderdg_kernel"].get("type","generic")=="optimised"
				useFortran   = useFortran or solver["aderdg_kernel"].get("language","C")=="Fortran"
			if "fv_kernel" in solver:
				useOptKernel = useOptKernel or solver["fv_kernel"].get("type","generic")=="optimised"
				useFortran   = useFortran or solver["fv_kernel"].get("language","C")=="Fortran"
		context["useOptKernel"] = useOptKernel
		context["useFortran"]   = useFortran
		
		generator.Generator.render_template(self,"Makefile.template",context,"Makefile","Generated ",True)
