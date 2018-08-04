from helper import BadSpecificationFile

class DirectoryAndPathChecker:
	def check(self, human_readable, pathInstance, required_subdirs=[], verbose_show_path=True, makedirs=False):
		valid = pathInstance.is_dir()
		human_status = "ok" if valid else "not found"
		if verbose_show_path or not valid:
			human_readable +=  ": " + pathInstance.resolve()
		if makedirs:
			if valid:
				human_status = "does exist (will not be overwritten)"
			else:
				
				os.makedirs(pathInstance) # todo, catch exceptions
				human_status = "created"

		if self.verbose:
			print(human_readable, " ... ", human_status)
		
		if not pathInstance.is_dir(): # check again in case of makedirs
			raise BadSpecificationFile("%s (relative to directory %s)" % (human_readable, os.getcwd()))
		
		for subdir in required_subdirs:
			check("%s holds %s",  pathInstance.joinpath(subdir), verbose_show_path=False, makedirs=makedirs)

	def __init__(self, spec, verbose):
		self.verbose = verbose
		
		paths = spec["paths"]
		peanoToolboxPath = paths["peanoKernelPath"] # sic, from old toolkit
		
		check("Peano kernel path", Path(paths["peanoKernelPath"]), required_subdirs=["tarch", "peano"])
		#check("Peano kernel path tarch sources", Path(paths["peanoKernelPath"]).joinpath("tarch/"))
		check("Peano toolboxes path", Path(peanoToolboxPath), required_subdirs=["multiscalelinkedcell", "sharedmemoryoracles", "mpibalancing"])
		
		check("ExaHyPE path", Path(paths["exahypePath"]))
		check("Output directory", Path(paths["outputDirectory"]), makedirs=True)
	
		# TODO (from Toolkit1): Initialize the CodeGeneratorHelper path (static)
		# CodeGeneratorHelper.setPaths(this);
		
