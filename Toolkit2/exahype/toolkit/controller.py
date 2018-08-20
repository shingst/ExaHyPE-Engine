#!/usr/bin/env python3

"""
This mimics the classical Frontend of the ExaHyPE toolkit.
"""

import os, sys, argparse, subprocess, datetime, json, logging
from pathlib import Path

# specfile module import
from specfiles import validate, OmniReader

# local import
from .directories import DirectoryAndPathChecker
from .toolkitHelper import BadSpecificationFile
from .toolkitHelper import ToolkitHelper
from .configuration import Configuration
from .solverController import SolverController
from .models import *


class Controller:
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
    
    def __init__(self):
        """TODO"""
        # Check the python version according to the configuration
        Configuration.checkPythonVersion()
        
        logging.basicConfig(format="%(filename)s:%(lineno)s(%(funcName)s):%(levelname)s %(message)s")
        self.log = logging.getLogger()

        # parse command line arguments
        args = self.parseArgs()
        
        # set member values from args
        self.specfileName = args.specfile.name # note, this can be something like "<stdin>"
        self.interactive = args.interactive or not args.not_interactive
        self.verbose = args.verbose or self.interactive
        self.write_json = args.write_json
        self.debug=args.debug
        
        if self.verbose:
            self.log.setLevel(logging.INFO)
            self.log.info(self.header())
            logging.raiseExceptions = True # development mode
            codegeneratorModel.CodegeneratorModel.logger = self.log #codegeneratorModel will print the associated command line
        else:
            logging.raiseExceptions = False # production mode
        
        inpath = Path(self.specfileName)
        if inpath.exists():
            self.log.info("Read input file %s." % inpath.resolve())
        else:
            self.log.info("Read from stdin (%s)" % str(args.specfile))
        
        rawSpec = self.getSpec(args.specfile, args.format)
        self.spec = self.validateAndSetDefaults(rawSpec, args.validate_only)


    def run(self):
        try:
            d = DirectoryAndPathChecker(self.buildBaseContext(), self.log)
        except BadSpecificationFile as e:
            self.log.error("Some directories did not exist and/or could not be created.")
            self.log.exception(e)
            self.log.error("Failure due to bad specificaiton file, cannot continue")
            sys.exit(-4)
        
        self.wait_interactive("validated and configured pathes")
        
        try:
            solverControl = SolverController(self.spec.get("solvers",[]), self.buildBaseContext())
            solverControl.run(self.log)
        except BadSpecificationFile as e:
            self.log.error("Could not create applications solver classes: %s" % str(e))
            self.log.exception(e)
            sys.exit(-6)
        
        self.wait_interactive("generated application-specific solver and plotter classes")
        
        try:
            # kernel calls
            kernelCalls = kernelCallsModel.KernelCallsModel(self.buildKernelCallsContext(codegeneratorModel.CodegeneratorModel.codegeneratorContextsList))
            pathToKernelCalls = kernelCalls.generateCode()
            self.log.info("Generated '"+pathToKernelCalls+"'")
        except Exception as e:
            self.log.error("Could not create ExaHyPE's kernel calls: %s" % str(e))
            self.log.exception(e)
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
            self.log.info("Generated '"+pathToMakefile+"'")
        except Exception as e:
            self.log.error("Could not create application-specific Makefile: %s" % str(e))
            self.log.exception(e)
            sys.exit(-10)
        
        self.wait_interactive("generated application-specific Makefile")
        
        self.checkEnvVariable()
        self.log.info(makefileMessage)
    
    
    def parseArgs(self):
        parser = argparse.ArgumentParser(
            description="ExaHyPE's new python-based toolkit",
            epilog="See http://www.exahype.eu and the Guidebook for more help."
        )
        
        # some mandatory arguments from Toolkit1
        optimized = parser.add_argument_group('Optimized kernels-specific options')
        optimized.add_argument("-c", "--clean-opt",
            help="Clean optimized kernels (only applicable if optimized kernels are used in Specfile)")

        ui = parser.add_argument_group("Toolkit user interface-specific options")
        g = ui.add_mutually_exclusive_group()
        g.add_argument("-i", "--interactive", action="store_true", default=False, help="Run interactively")
        g.add_argument("-n", "--not-interactive", action="store_true", default=True, help="Run non-interactive in non-stop-mode (default)")
        ui.add_argument("-v", "--verbose", action="store_true", help="Be verbose (off by default, triggered on by --interactive)")
        ui.add_argument("-j", "--write-json", action="store_true", default=False, help="Write a JSON file if legacy specification file is read (*.exahype) (default: no)")
        ui.add_argument("-d", "--debug", action="store_true", default=False, help="Turn the debug mode on: print stack traces and more detailed error messages.")

        formats = parser.add_argument_group("Specification-file-format specific options")
        formats.add_argument("-o", "--validate-only", action="store_true", default=False, help="Validate input only, don't run the toolkit. Will output the correct JSON if passes.")
        formats.add_argument("-f", "--format", choices=OmniReader.available_readers(), default=OmniReader.any_format_name, help="Specification file format of the input file. 'any' will try all built-in-formats.")

        parser.add_argument('specfile',
            type=argparse.FileType('r'),
            help="The specification file to work on (can be .exahype, .exahype2, .json, etc.)")
        
        return parser.parse_args()
    
    
    def getSpec(self, file_handle, file_format=OmniReader.any_format_name):
        try:
            reader = OmniReader(self.log)
            spec = reader.read(file_handle.read(), required_file_format=file_format, filename=self.specfileName)
        except Exception as e:
            self.log.error("Could not read specification file '%s': %s" % (self.specfileName, str(e)))
            self.log.error("In order to fix this problem, please fix the format of your file with the command line flag --format=XXX where XXX is a supported specification file format.")
            if self.debug:
               self.log.exception(e)
            sys.exit(-3)

        # I find this is more a debugging feature...
        if self.specfileName.endswith(".exahype") and self.write_json:
            json_file_name = self.specfileName.replace(".exahype",".exahype2")
            with open(json_file_name, 'w') as outfile:
                json.dump(spec,outfile,indent=2)
                self.log.info("Write JSON file '%s' ... OK" % json_file_name)

        return spec
    
    def validateAndSetDefaults(self, spec, validate_only=False):
        """
        Given a specification, validate it  against the JSON-Schema and
        returns the native python data structure, enriched with default values from the schema.
        """
        try:
            spec = validate(spec, set_defaults=True)
        except Exception as e:
            self.log.error("Specification file does not hold a valid ExaHyPE specification, it did not pass the schema validation step. The error message is: %s" % str(e))
            self.log.exception(e)
            sys.exit(-4)

        if validate_only:
            print(json.dumps(spec, sort_keys=True, indent=4))
            sys.exit(0)
        else:
            return spec

    def buildBaseContext(self):
        """Generate base context from spec with commonly used value"""
        context = {}
        # commonly used paths
        context["outputPath"]           = os.path.join(Configuration.pathToExaHyPERoot, self.spec["paths"]["output_directory"])
        context["exahypePath"]          = os.path.join(Configuration.pathToExaHyPERoot, self.spec["paths"]["exahype_path"])
        context["peanoToolboxPath"]     = os.path.join(Configuration.pathToExaHyPERoot, self.spec["paths"]["peano_kernel_path"])
        context["plotterSubDirectory"]  = self.spec["paths"].get("plotter_subdirectory","").strip()
        
        # commonly used parameters
        context["project"]          = self.spec["project_name"]
        
        context["architecture"]     = self.spec["architecture"]
        context["alignment"]        = Configuration.alignmentPerArchitectures[context["architecture"]]
        
        context["dimensions"]   = self.spec["computational_domain"]["dimension"]
        context["range_0_nDim"] = range(0,context["dimensions"])
        
        context["enableProfiler"] = False # TODO read from spec
        
        return context
    

    def buildMakefileContext(self):
        """Generate context for the Makefile model"""
        context = self.buildBaseContext()
        
        context["architecture"]      = self.spec["architecture"]
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
    
    def buildKernelCallsContext(self, codegeneratorContextsList):
        """Generate context for the KernelCalls model"""
        context = self.buildBaseContext()
        
        context["solvers"] = []
        for solver in self.spec.get("solvers",[]):
            solverContext = {}
            solverContext["name"]                        = solver["name"]
            solverContext["type"]                        = solver["type"]
            solverContext["class"]                       = context["project"]+"::"+solver["name"]
            solverContext["headerPath"]                  = solver["name"]+".h"
            solverContext["variables_as_str"]            = ToolkitHelper.variables_to_str(solver,"variables")
            solverContext["material_parameters_as_str"]  = ToolkitHelper.variables_to_str(solver,"material_parameters")
            solverContext["global_observables_as_str"]   = ToolkitHelper.variables_to_str(solver,"global_observables")
            solverContext["plotters"] = []
            for plotter in solver.get("plotters",[]):
                plotterContext = {}
                plotterContext["name"]             = plotter["name"]
                plotterContext["headerPath"]       = os.path.join(context["plotterSubDirectory"], plotter["name"]+".h" )
                plotterContext["type_as_str"]      = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
                plotterContext["variables_as_str"] = ToolkitHelper.variables_to_str(plotter,"variables")
                solverContext["plotters"].append(plotterContext)
            context["solvers"].append(solverContext)
        
        context["specfileName"]          = self.specfileName
        context["specFileAsHex"]         = self.specfileAsHex(self.spec)
        context["externalParserCommand"] = "%s/%s %s" % ( Configuration.pathToExaHyPERoot, "Toolkit2/toolkit.sh","--format=any --validate-only")
        context["codegeneratorContextsList"] = codegeneratorContextsList #list of contexts used for the generation of each optimized solver
        # todo(JM) profiler 
        # todo(Sven) serialised spec file compiled into KernelCalls.cp
        context["includePaths"] = [] #TODO
        
        return context
    
    def specfileAsHex(self,spec):
        """
        Given a native python nested dict/list object, dump it as string and then hex-encode that string
        character by character. This is safest way to include something in C++ without dealing with
        character sets or anything.
        """
        text = json.dumps(spec, sort_keys=True, indent=4)
        hex_tokens = [ "0x%02x"%ord(char) for char in text ] + ["0x00"] # null-terminated list of hex numbers
        return ", ".join(hex_tokens)
    
    def checkEnvVariable(self):
        """
        Check environment variables as an ultimate step. Should only be called when actually calling the toolkit,
        so it is distinct from the validation step.
        """
        if "shared_memory" in self.spec:
            if not "TBB_INC" in os.environ:
                self.log.warning("environment variable TBB_INC not set but required if code is built with TBB");
            if not "TBB_SHLIB" in os.environ:
                self.log.warning("environment variable TBB_SHLIB not set but required if code is built with TBB");
