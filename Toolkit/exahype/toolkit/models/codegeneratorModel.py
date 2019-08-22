
import sys

# add path to dependencies
from ..configuration import Configuration
sys.path.insert(1, Configuration.pathToCodegenerator)
import codegenerator


class CodegeneratorModel:


    def generateCode(self, solverContext):
        #translation from solverContext to codegeneratorContext
        if solverContext["kernelType"] == "aderdg":
            codegeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "aderdg",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "order"              : solverContext["order"],
                "dimension"          : solverContext["dimensions"],
                "numerics"           : "linear" if solverContext["isLinear"] else "nonlinear",
                "architecture"       : solverContext["architecture"],
                # Optional bool parameters (may set redundant flags and default false flag)
                "useFlux"            : solverContext["useFlux"],
                "useFluxVect"        : solverContext["useFluxVect"],
                "useNCP"             : solverContext["useNCP"],
                "useNCPVect"         : solverContext["useNCPVect"],
                "useSource"          : solverContext["useSource"],
                "useSourceVect"      : solverContext["useSourceVect"],
                "useFusedSource"     : solverContext["useFusedSource"],
                "useFusedSourceVect" : solverContext["useFusedSourceVect"],
                "useMaterialParam"   : solverContext["useMaterialParameters"],
                "useMaterialParamVect" : solverContext["useMaterialParametersVect"],
                "useCERKGuess"       : solverContext["useCERK"],
                "useSplitCKScalar"   : solverContext["useSplitCKScalar"],
                "useSplitCKVect"     : solverContext["useSplitCKVect"],
                "useGaussLobatto"    : solverContext["basis"] == "lobatto",
                # Optional int parameters (may set redundant flags)
                "usePointSources"    : solverContext["numberOfPointSources"] if solverContext["numberOfPointSources"] > 0 else -1,
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        elif solverContext["kernelType"] == "limiter":
            codegeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "limiter",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "order"              : solverContext["order"],
                "dimension"          : solverContext["dimensions"],
                "architecture"       : solverContext["architecture"],
                "numberOfObservable" : solverContext.get("numberOfDMPObservables", -1), #not set if not limiterSolver
                "ghostLayerWidth"    : solverContext.get("ghostLayerWidth", 0), #not set if not limiterSolver
                "limPatchSize"       : solverContext.get("patchSize",-1),
                # Optional bool parameters (may set redundant flags and default false flag)
                "useGaussLobatto"    : solverContext["basis"] == "lobatto",
                # Optional int parameters (may set redundant flags)
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        elif solverContext["kernelType"] == "fv":
            codegeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "fv",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "patchSize"          : solverContext["patchSize"],
                "dimension"          : solverContext["dimensions"],
                "finiteVolumesType"  : solverContext["finiteVolumesType"],
                "architecture"       : solverContext["architecture"],
                # Optional bool parameters (may set redundant flags and default false flag)
                "useFlux"            : solverContext["useFlux"],
                "useViscousFlux"     : solverContext["useViscousFlux"],
                "useNCP"             : solverContext["useNCP"],
                "useSource"          : solverContext["useSource"],
                "useFusedSource"     : solverContext["useFusedSource"],
                "useMaterialParam"   : solverContext["useMaterialParameters"],
                # Optional int parameters (may set redundant flags)
                "usePointSources"    : solverContext["numberOfPointSources"] if solverContext["numberOfPointSources"] > 0 else -1,
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        else:
            raise ValueError("KernelType '"+context["kernelType"]+"' not supported by the codegenerator")
        # call the codegenerator with the given context
        codegeneratorController = codegenerator.Controller(codegeneratorContext)
        codegeneratorController.generateCode()
        
        # if verbose print the associated command line
        codegeneratorContext["commandLine"] = codegeneratorController.commandLine
        
        return solverContext["optKernelPath"], codegeneratorContext
