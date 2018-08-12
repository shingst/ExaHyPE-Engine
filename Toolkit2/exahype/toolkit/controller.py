#!/usr/bin/env python3

"""
This mimics the classical Frontend of the ExaHyPE toolkit.
"""

import os, sys, argparse, subprocess, datetime, json

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
        self.write_json=args.write_json
        self.debug_mode=args.debug
        
        if self.verbose: # otherwise no need to call git, etc. 
            self.info(self.header())
        
        inpath = Path(self.specfileName)
        if inpath.exists():
            self.info("Read input file %s." % inpath.resolve())
        else:    self.info("Read from stdin (%s)" % str(args.specfile))
        
        self.spec = self.getSpec(args.specfile)
    
    
    def run(self):
        try:
            d = directories.DirectoryAndPathChecker(self.buildBaseContext(), self.verbose)
        except BadSpecificationFile as e:
            print("ERROR: Some directories did not exist and/or could not be created.", file=sys.stderr)
            print("ERROR: Reason: %s" % e, file=sys.stderr)
            print("ERROR: Failure due to bad specificaiton file, cannot continue", file=sys.stderr)
            print("\n\n** Completed with errors **")
            sys.exit(-4)
        
        self.wait_interactive("validated and configured pathes")
        
        try:
            for i,solver in enumerate(self.spec.get("solvers",[])):
                print("Generating solver[%d] = %s..." % (i, solver["name"]))
                solverFiles = []
                if solver["type"]=="ADER-DG":
                    model       = solverModel.SolverModel(self.buildADERDGSolverContext(solver))
                    solverFiles = model.generateCode()
                elif solver["type"]=="Finite-Volumes":
                    model       = solverModel.SolverModel(buildFVSolverContext(solver))
                    solverFiles = model.generateCode()
                elif solver["type"]=="Limiting-ADER-DG":
                    aderdgContext = self.buildADERDGSolverContext(solver)
                    fvContext     = self.buildFVSolverContext(solver)
                    context       = self.buildLimitingADERDGSolverContext(solver)
                    # modifications
                    fvContext["solver"]         = context["FVSolver"]
                    fvContext["abstractSolver"] = context["AbstractFVSolver"]
                    fvContext["patchSize"]      = 2 * aderdgContext["order"] + 1
                    
                    aderdgContext["solver"]                 = context["ADERDGSolver"]
                    aderdgContext["abstractSolver"]         = context["AbstractADERDGSolver"]
                    aderdgContext["numberOfDMPObservables"] = context["numberOfDMPObservables"]
                    
                    # generate all solver files
                    model        = solverModel.SolverModel(context)
                    solverFiles  = model.generateCode()
                    model        = solverModel.SolverModel(fvContext)
                    solverFiles += model.generateCode()
                    model        = solverModel.SolverModel(aderdgContext)
                    solverFiles += model.generateCode()
                
                for path in solverFiles:
                    if not path is None:
                        print("Generated '"+path+"'")
                for j,plotter in enumerate(solver.get("plotters",[])):
                    print("Generating plotter[%d] = %s for solver[%d] = %s" % (j, plotter["name"], i, solver["name"]))
                    model = plotterModel.PlotterModel(self.buildPlotterContext(solver,plotter))
                    for path in model.generateCode():
                        if not path is None:  
                            print("Generated '"+path+"'")
        except BadSpecificationFile as e:
            print("ERROR: Could not create applications solver classes.", file=sys.stderr)
            print("ERROR: Reason: %s" % e,file=sys.stderr)
            print("\n\n** Completed with errors **")
            sys.exit(-6)
        
        self.wait_interactive("generated application-specific solver and plotter classes")
        
        try:
            # kernel calls
            kernelCalls = kernelCallsModel.KernelCallsModel(self.buildKernelCallsContext())
            pathToKernelCalls = kernelCalls.generateCode()
            print("Generated '"+pathToKernelCalls+"'")
        except Exception as e:
            print("ERROR: Could not create ExaHyPE's kernel calls", file=sys.stderr)
            print("ERROR: Reason: %s" % e,file=sys.stderr)
            print("\n\n** Completed with errors **")
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
            print("Generated '"+pathToMakefile+"'")
        except Exception as e:
            print("ERROR Could not create application-specific Makefile", file=sys.stderr)
            print("ERROR: Reason: %s" % e,file=sys.stderr)
            print("\n\n** Completed with errors **")
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
        parser.add_argument("-j", "--write-json", action="store_true", default=False, help="Write a JSON file if legacy specification file is read (*.exahype) (default: no)")
        parser.add_argument("-d", "--debug", action="store_true", default=False, help="Turn the debug mode on: print stack traces and more detailed error messages.")
        parser.add_argument('specfile',
            type=argparse.FileType('r'),
            help="The specification file to work on (can be .exahype, .exahype2, .json)")
        
        return parser.parse_args()
    
    
    def getSpec(self, specfilePath):
        try:
            return self.load(self.specfileName)
        except Exception as e:
            print("ERROR: Could not properly read specfile",file=sys.stderr)
            print("ERROR: Reason: %s" % e,file=sys.stderr)
            print("\n\n** Completed with errors **")
            if self.debug:
              raise
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
          if self.write_json:
            json_file_name = specfile_name.replace(".exahype",".exahype2")
            with open(json_file_name, 'w') as outfile:
              json.dump(specification,outfile,indent=2)
              print("Write JSON file '%s' ... OK" % json_file_name )
        else:
          specification = validate(specfile_name)
        return specification
    
    
    def buildBaseContext(self):
        """Generate base context from spec with commonly used value"""
        context = {}
        # commonly used paths
        context["outputPath"]           = Configuration.pathToExaHyPERoot+"/"+self.spec["paths"]["output_directory"]
        context["exahypePath"]          = Configuration.pathToExaHyPERoot+"/"+self.spec["paths"]["exahype_path"]
        context["peanoToolboxPath"]     = Configuration.pathToExaHyPERoot+"/"+self.spec["paths"]["peano_kernel_path"]
        context["plotter_subdirectory"] = self.spec["paths"]["peano_kernel_path"]
        
        # commonly used parameters
        context["project"]          = self.spec["project_name"]
        
        context["alignment"]        = Configuration.alignmentPerArchitectures[self.spec["architecture"]]
        
        context["dimensions"]   = self.spec["computational_domain"]["dimension"]
        context["range_0_nDim"] = range(0,context["dimensions"])
        
        context["enableProfiler"] = False # TODO read from spec
        
        return context
    
    def buildADERDGSolverContext(self,solver):
        context = self.buildBaseContext()
        context.update(self.buildBaseSolverContext(solver))
        context.update(self.buildADERDGKernelContext(solver["aderdg_kernel"]))
        
        context["order"]                  = solver["order"]
        context["numberOfDMPObservables"] = 0 # overwrite if called from LimitingADERDGSolver creation
        
        return context
        
    def buildLimitingADERDGSolverContext(self,solver):
        context = self.buildBaseContext()
        context.update(self.buildBaseSolverContext(solver))
        
        context["order"]                  = solver["order"]
        context["numberOfDMPObservables"] = solver["limiter"]["numberOfObservables"]
        context["implementation"]         = solver["limiter"].get("implementation","generic")
        
        context["ADERDGSolver"]         = solver["name"]+"_ADERDG"
        context["FVSolver"]             = solver["name"]+"_FV"
        context["ADERDGAbstractSolver"] = "Abstract"+solver["name"]+"_ADERDG"
        context["FVAbstractSolver"]     = "Abstract"+solver["name"]+"_FV"
        
        return context
    
        
    def buildFVSolverContext(self,solver):
        context = self.buildBaseContext()
        context.update(self.buildBaseSolverContext(solver))
        context.update(self.buildFVKernelContext(solver["fv_kernel"]))
        
        context["patchSize"] = solver["patch_size"] # overwrite if called from LimitingADERDGSolver creation
        
        return context
    
    def buildBaseSolverContext(self,solver):
        context = {}
        
        context["solverType"]     = solver["type"]
        context["solver"]         = solver["name"]
        context["abstractSolver"] = "Abstract"+context["solver"]
        
        nVar          = helper.count_variables(helper.parse_variables(solver,"variables"))
        nParam        = helper.count_variables(helper.parse_variables(solver,"material_parameters"))
        nGlobalObs    = helper.count_variables(helper.parse_variables(solver,"global_observables"))
        nPointSources = len(solver.get("point_sources",[]))
        
        context["numberOfVariables"]          = nVar
        context["numberOfMaterialParameters"] = nParam
        context["numberOfGlobalObservables"]  = nGlobalObs
        context["numberOfPointSources"]       = nPointSources
        
        context["range_0_nVar"]          = range(0,nVar)
        context["range_0_nVarParam"]     = range(0,nVar+nParam)
        context["range_0_nGlobalObs"]    = range(0,nGlobalObs)    # nGlobalObs might be 0
        context["range_0_nPointSources"] = range(0,nPointSources) # nPointSources might be 0
        
        context["namingSchemes"]=[] # TODO read from spec
        return context
    
    def buildADERDGKernelContext(self,kernel):
        context = {}
        context["implementation"]          = kernel.get("implementation","generic")
        context["useMaxPicardIterations"]  = kernel.get("space_time_predictor",{}).get("maxpicarditer",0)!=0 
        context["maxPicardIterations"]     = kernel.get("space_time_predictor",{}).get("maxpicarditer",0)
        context["tempVarsOnStack"]         = kernel.get("allocate_temporary_arrays","heap")=="stack" 
        context["patchwiseAdjust"]         = kernel.get("adjust_solution","pointwise")=="patchwise" 
        context["language"]                = kernel.get("language","C").lower()
        context["basis"]                   = kernel.get("basis","Legendre").lower()
        context["isLinear"]                = not kernel.get("nonlinear",True)
        context["isNonlinear"]             = kernel.get("nonlinear",True)
        context["linearOrNonlinear"]       = "Linear" if context["isLinear"] else "Nonlinear"
        context["isFortran"]               = kernel.get("language",False)=="Fortran" 
        context["useCERK"]                 = kernel.get("space_time_predictor",{}).get("cerkguess",False)
        context["noTimeAveraging"]         = "true" if kernel.get("space_time_predictor",{}).get("notimeavg",False) else "false"
        context["useConverter"]            = "converter" in kernel.get("optimised_kernel_debugging",[])
        context["countFlops"]              = "flops" in kernel.get("optimised_kernel_debugging",[])
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context
        
    def buildFVKernelContext(self,kernel):
        context = {}
        context["implementation"] = kernel.get("implementation","generic")
        ghostLayerWidth = { "godunov" : 1, "musclhancock" : 2 }
        context["ghostLayerWidth"]=ghostLayerWidth[kernel["scheme"]]
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context
    
    def buildKernelTermsContext(self,terms):
        context = {}
        for term in ["flux","source","ncp","point_sources","material_parameters"]:
            option = term.replace("_s","S").replace("_p","P").replace("ncp","NCP")
            option = "use%s%s" % ( option[0].upper(), option[1:] )
            context[option] = term in terms
            context["%s_s" % option] = "true" if context[option] else "false"
        return context
        
    def buildPlotterContext(self,solver,plotter):
        context = self.buildBaseContext()
        context.update(self.buildBaseSolverContext(solver))

        context["plotter"]         = plotter["name"]
        context["writtenUnknowns"] = helper.count_variables(helper.parse_variables(solver,"variables"))
        context["plotterType"]     = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
        
        context["plotterDir"]=""
        subdirOption = "plotter_subdirectory"
        if subdirOption in self.spec["paths"] and len(self.spec["paths"][subdirOption].strip()):
            context["plotterDir"]  = self.spec["paths"][subdirOption]+"/"
            context["outputPath"] += self.spec["paths"][subdirOption]
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
    
    def buildKernelCallsContext(self):
        """Generate context for the KernelCalls model"""
        context = self.buildBaseContext()
        
        context["solvers"] = []
        plotter_subdirectory = self.spec["paths"].get("plotter_subdirectory","").strip()
        for solver in self.spec.get("solvers",[]):
            solverContext = {}
            solverContext["name"]                        = solver["name"]
            solverContext["type"]                        = solver["type"]
            solverContext["class"]                       = context["project"]+"::"+solver["name"]
            solverContext["headerPath"]                  = solver["name"]+".h"
            solverContext["variables_as_str"]            = helper.variables_to_str(solver,"variables")
            solverContext["material_parameters_as_str"]  = helper.variables_to_str(solver,"material_parameters")
            solverContext["global_observables_as_str"]   = helper.variables_to_str(solver,"global_observables")
            solverContext["plotters"] = []
            for plotter in solver.get("plotters",[]):
                plotterContext = {}
                plotterContext["name"]             = plotter["name"]
                plotterContext["headerPath"]       = os.path.join(plotter_subdirectory, plotter["name"]+".h" )
                plotterContext["type_as_str"]      = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
                plotterContext["variables_as_str"] = helper.variables_to_str(plotter,"variables")
                solverContext["plotters"].append(plotterContext)
            context["solvers"].append(solverContext)
        
        context["specfileName"]     = self.specfileName
        context["spec_file_as_hex"] = "0x2F" # todo(Sven) do the conversion
        context["subPaths"]         = []
        # todo(JM) optimised kernels subPaths
        # todo(JM) profiler 
        # todo(Sven) serialised spec file compiled into KernelCalls.cp
        context["includePaths"] = [] #TODO
        
        return context
    
    #TODO shouldn't this be in the validation step?
    def checkEnvVariable(self):
        if "shared_memory" in self.spec:
            if not "TBB_INC" in os.environ:
                print("WARNING: environment variable TBB_INC not set but required if code is built with TBB");
            if not "TBB_SHLIB" in os.environ:
                print("WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB");
