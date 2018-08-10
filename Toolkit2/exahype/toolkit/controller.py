#!/usr/bin/env python3

"""
This mimics the classical Frontend of the ExaHyPE toolkit.
"""

import os, sys, argparse, subprocess, datetime

from os.path import isdir, isfile
from pathlib import Path
from collections import OrderedDict

sys.path.append(os.path.join(os.path.dirname(__file__),"..","..")) # to allow import exahype... work
from exahype.toolkit import *
from exahype.toolkit.helper import BadSpecificationFile
from exahype.specfiles import validate,validate_specfile1

from models import *
from configuration import Configuration

class Controller():
    def header(self):
        info = {
            "gittag": subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip().decode('ascii'),
            "year": str(datetime.datetime.now().year)
        }
        return """
    ______           __  __      ____  ______    
   / ____/  ______ _/ / / /_  __/ __ \/ ____/    *************************
  / __/ | |/_/ __ `/ /_/ / / / / /_/ / __/        The ExaHyPE Toolkit v2
 / /____>  </ /_/ / __  / /_/ / ____/ /___           www.exahype.eu
/_____/_/|_|\__,_/_/ /_/\__, /_/   /_____/          Commit: %(gittag)s
                       /____/                    *************************

 This project has received funding from the European Union's Horizon 2020
 research and innovation programme under grant agreement No 671698 (ExaHyPE).

 Copyright (c) %(year)s, All rights reserved
 ExaHyPE is based  on the PDE framework Peano (www.peano-framework.org).
 Released under the BSD 3 Open Source License.
 
""" % info
        # end of header

    def wait_interactive(self, msg):
        if self.interactive:
            print("\n\n\n\n")
            print(msg + "... ok")
            input("<pres Enter>")
            print("\n\n\n\n")
            
    def info(self, msg):
        if self.verbose:
            print(msg)
    
    
    def __init__(self):
        # parse command line arguments
        args = self.parseArgs()
        
        #set member values from args
        self.specfileName = args.specfile.name
        self.interactive = args.interactive or not args.not_interactive
        self.verbose = args.verbose or self.interactive
        
        if self.verbose: # otherwise no need to call git, etc. 
            self.info(self.header())
        
        inpath = Path(self.specfileName)
        if inpath.exists():
            self.info("Read input file %s." % inpath.resolve())
        else:    self.info("Read from stdin (%s)" % str(args.specfile))
        
        self.spec = self.getSpec(args.specfile)
    
    
    def run(self):
        try:
            d = directories.DirectoryAndPathChecker(self.spec, self.verbose)
        except BadSpecificationFile as e:
            print("Some directories did not exist and/or could not be created.")
            print("Error message: ", e)
            print("Failure due to bad specificaiton file, cannot continue")
            sys.exit(-4)
        
        self.wait_interactive("validated and configured pathes")
        
        try:
            s = solver.SolverGenerator(self.spec, self.specfileName, self.verbose)
            s.generate_all_solvers()
        except BadSpecificationFile as e:
            print("Could not create applications solver classes.")
            print(e)
            sys.exit(-6)
        
        self.wait_interactive("generated application-specific solver classes")
        
        try:
            # kernel calls
            kernelCalls = kernelCallsModel.KernelCallsModel(self.buildKernelCallsContext())
            pathToKernelCalls = kernelCalls.generateCode()
            print("Generated "+pathToKernelCalls)
        except Exception as e:
            print("Could not create ExaHyPE's kernel calls")
            print(e)
            sys.exit(-10)
            
        self.wait_interactive("generated computational kernel calls")
        
        # blah
        #CreateReadme(self.spec,self.verbose)
        #
        #self.wait_interactive("generated README")
        #
        # makefiles, etc.
        #setupBuildEnvironment(self.spec, self.verbose)
        #
        #self.wait_interactive("setup build environment")
        
        try:
            makefile = makefileModel.MakefileModel(self.buildMakefileContext())
            pathToMakefile = makefile.generateCode() #generate Makefile and get path to it
            makefileMessage = makefile.getOutputMessage()
            print("Generated "+pathToMakefile)
        except Exception as e:
            print("Could not create application-specific Makefile")
            print(e)
            sys.exit(-10)
        
        self.wait_interactive("generated application-specific Makefile")
        
        self.checkEnvVariable()
        print(makefileMessage)
    
    
    def parseArgs(self):
        parser = argparse.ArgumentParser(
            description="ExaHyPE's new python-based toolkit",
            epilog="See http://www.exahype.eu and the Guidebook for more help."
        )
        
        # some mandatory arguments from Toolkit1
        parser.add_argument("-c", "--clean-opt",
            help="Clean optimized kernels (only applicable if optimized kernels are used in Specfile)")
        g = parser.add_mutually_exclusive_group()
        g.add_argument("-i", "--interactive", action="store_true", default=False, help="Run interactively")
        g.add_argument("-n", "--not-interactive", action="store_true", default=True, help="Run non-interactive in non-stop-mode (default)")
        parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose (off by default, triggered on by --interactive)")

        parser.add_argument('specfile',
            type=argparse.FileType('r'),
            help="The specification file to work on (can be .exahype, .exahype2, .json)")
        
        return parser.parse_args()
    
    
    def getSpec(self, specfilePath):
        try:
            return self.load(self.specfileName)
        except Exception as e:
            print("Could not properly read specfile")
            print(e)
            sys.exit(-3)
    
    
    def load(self, specfile_name):
        """
        Given a specfile name, it assumes it to be JSON, reads it, validates it
        against the JSON-Schema and returns the native python data structure,
        enriched with default values from the schema.
        """
        specification=""
        if specfile_name.endswith(".exahype"):
          specification = validate_specfile1(specfile_name)
        else:
          specification = validate(specfile_name)
        return specification
    
    
    def buildBaseContext(self):
        """Generate base context from spec with commonly used value"""
        context = {}
        # outputPath for generated file
        context["outputPath"] = self.spec["paths"]["output_directory"]
        # commonly used values
        context["dimensions"] = self.spec["computational_domain"]["dimension"]
        context["project"] = self.spec["project_name"]
        context["exahypePath"] = self.spec["paths"]["exahype_path"]
        context["peanoToolboxPath"] = self.spec["paths"]["peano_kernel_path"]
        context["alignment"] = Configuration.alignmentPerArchitectures[self.spec["architecture"]]
        
        return context
    
    
    def buildMakefileContext(self):
        """Generate context for the Makefile model"""
        context = self.buildBaseContext()
        context["useSharedMem"]      = "shared_memory" in self.spec;
        context["useDistributedMem"] = "distributed_memory" in self.spec;
        context["useIpcm"]   = False # TODO
        context["useLikwid"] = False # TODO
        context["likwidInc"] = ""    # TODO
        # kernels
        useOptKernel = False
        useFortran   = False
        for solver in self.spec["solvers"]:
            if "aderdg_kernel" in solver:
                useOptKernel = useOptKernel or solver["aderdg_kernel"].get("type","generic")=="optimised"
                useFortran   = useFortran or solver["aderdg_kernel"].get("language","C")=="Fortran"
            if "fv_kernel" in solver:
                useOptKernel = useOptKernel or solver["fv_kernel"].get("type","generic")=="optimised"
                useFortran   = useFortran or solver["fv_kernel"].get("language","C")=="Fortran"
        context["useOptKernel"] = useOptKernel
        context["useFortran"]   = useFortran
        
        return context
    
    def buildKernelCallsContext(self):
        """Generate context for the KernelCalls model"""
        context = self.buildBaseContext()
        
        context["solvers"] = []
        plotter_subdirectory = self.spec["paths"].get("plotter_subdirectory","").strip()
        for solver in self.spec.get("solvers",[]):
            solverContext = {}
            solverContext["name"] = solver["name"]
            solverContext["type"] = solver["type"]
            solverContext["class"] = context["project"]+"::"+solver["name"]
            solverContext["headerPath"] = solver["name"]+".h"
            solverContext["variables_as_str"]            = helper.variables_to_str(solver,"variables")
            solverContext["material_parameters_as_str"]  = helper.variables_to_str(solver,"material_parameters")
            solverContext["global_observables_as_str"]   = helper.variables_to_str(solver,"global_observables")
            solverContext["plotters"] = []
            for plotter in solver.get("plotters",[]):
                plotterContext = {}
                plotterContext["headerPath"]      = os.path.join(plotter_subdirectory, plotter["name"]+".h" )
                plotterContext["type_as_str"]      = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
                plotterContext["variables_as_str"] = helper.variables_to_str(plotter,"variables")
                solverContext["plotters"].append(plotterContext)
            context["solvers"].append(solverContext)
        
        context["specfileName"]        = self.specfileName
        context["spec_file_as_hex"] = "0x2F" # todo(Sven) do the conversion
        context["subPaths"]         = []
        # todo(JM) optimised kernels subPaths
        # todo(JM) profiler 
        # todo(Sven) serialised spec file compiled into KernelCalls.cp
        context["includePaths"] = []
        
        return context
    
    #TODO shouldn't this be in the validation step?
    def checkEnvVariable(self):
        if "shared_memory" in self.spec:
            if not "TBB_INC" in os.environ:
                print("WARNING: environment variable TBB_INC not set but required if code is built with TBB");
            if not "TBB_SHLIB" in os.environ:
                print("WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB");
