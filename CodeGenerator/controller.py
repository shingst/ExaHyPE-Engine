#!/bin/env python
##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# Controller of the code generator
#
# @note
# requires python3


import argparse
import os
import copy
import subprocess
import errno
import time

from models import *


class Controller:
    """Main Controller
    
    Read the input from the command line, validate them and generate a base 
    context for the models.
    
    Use generateCode() to run the models with the base context.
    
    Can generate gemms with generateGemms(outputFile, matmulconfig), will be done
    automatically when using generateCode().
    """
    
    def __init__(self, absolutePathToRoot, absolutePathToLibxsmm, simdWidth):
        """Initialize the base config from the command line inputs"""
        # Process the command line arguments
        parser = argparse.ArgumentParser(description="This is the front end of the ExaHyPE code generator.")
        # Mandatory parameters
        parser.add_argument("pathToApplication",
            help="path to the application as given by the ExaHyPE specification file (application directory as root)")
        parser.add_argument("pathToOptKernel",
            help="desired relative path to the generated code (application directory as root)")
        parser.add_argument("namespace",
            help="desired namespace for the generated code")
        parser.add_argument("solverName",
            help="name of the user-solver")
        parser.add_argument("numberOfVariables",
            type=int,
            help="the number of quantities")
        parser.add_argument("numberOfParameters",
            type=int,
            help="the number of parameters (fixed quantities)")
        parser.add_argument("order",
            type=int,
            help="the order of the approximation polynomial")
        parser.add_argument("dimension",
            type=int,
            help="number of dimensions you want to simulate")
        parser.add_argument("numerics",
            help="linear or nonlinear")
        parser.add_argument("architecture",
            help="the microarchitecture of the target device")
        # Optional parameters
        parser.add_argument("--useFlux",
            action="store_true",
            help="enable flux")
        parser.add_argument("--useFluxVect",
            action="store_true",
            help="enable vectorized flux (include useFlux)")
        parser.add_argument("--useNCP",
            action="store_true",
            help="enable non conservative product")
        parser.add_argument("--useNCPVect",
            action="store_true",
            help="enable vectorized non conservative product  (include useNCP)")
        parser.add_argument("--useSource",
            action="store_true",
            help="enable source terms")
        parser.add_argument("--useSourceVect",
            action="store_true",
            help="enable vectorized source terms (include useSource)")
        parser.add_argument("--useFusedSource",
            action="store_true",
            help="enable fused source terms (include useSource)")
        parser.add_argument("--useFusedSourceVect",
            action="store_true",
            help="enable vectorized fused source terms (include useFusedSource and useSourceVect)")
        parser.add_argument("--useMaterialParam",
            action="store_true",
            help="enable material parameters")
        parser.add_argument("--usePointSources",
            type=int,
            default=-1, #default -1 marks option not used
            metavar='nPointSources',
            help="enable nPointSources point sources")
        parser.add_argument("--useCERKGuess",
            action="store_true",
            help="use CERK for SpaceTimePredictor inital guess (nonlinear only)")
        parser.add_argument("--useLimiter",
            type=int,
            default=-1, #default -1 marks option not used
            metavar='useLimiter',
            help="enable limiter with the given number of observable")
        parser.add_argument("--useGaussLobatto",
            action="store_true",
            help="use Gauss Lobatto Quadrature instead of Gauss Legendre")
        parser.add_argument("--ghostLayerWidth",
            type=int,
            default=0,
            metavar='ghostLayerWidth',
            help="use limiter with the given ghostLayerWidth, requires useLimiter option, default = 0")
        commandLineArguments = parser.parse_args()
        
        # Generate the base config from the parsed input
        self.config = {
            "numerics"              : commandLineArguments.numerics,
            "pathToOptKernel"       : commandLineArguments.pathToOptKernel,
            "solverName"            : commandLineArguments.solverName,
            "nVar"                  : commandLineArguments.numberOfVariables,
            "nPar"                  : commandLineArguments.numberOfParameters,
            "nData"                 : commandLineArguments.numberOfVariables + commandLineArguments.numberOfParameters,
            "nDof"                  : (commandLineArguments.order)+1,
            "nDim"                  : commandLineArguments.dimension,
            "useFlux"               : (commandLineArguments.useFlux or commandLineArguments.useFluxVect),
            "useFluxVect"           : commandLineArguments.useFluxVect,
            "useNCP"                : (commandLineArguments.useNCP or commandLineArguments.useNCPVect),
            "useNCPVect"            : commandLineArguments.useNCPVect,
            "useSource"             : (commandLineArguments.useSource or commandLineArguments.useSourceVect or commandLineArguments.useFusedSource or commandLineArguments.useFusedSourceVect),
            "useSourceVect"         : commandLineArguments.useSourceVect,
            "useFusedSource"        : (commandLineArguments.useFusedSource or commandLineArguments.useFusedSourceVect),
            "useFusedSourceVect"    : commandLineArguments.useFusedSourceVect,
            "nPointSources"         : commandLineArguments.usePointSources,
            "usePointSources"       : commandLineArguments.usePointSources >= 0,
            "useMaterialParam"      : commandLineArguments.useMaterialParam,
            "codeNamespace"         : commandLineArguments.namespace,
            "pathToOutputDirectory" : os.path.join(absolutePathToRoot,commandLineArguments.pathToApplication,commandLineArguments.pathToOptKernel),
            "architecture"          : commandLineArguments.architecture,
            "useLimiter"            : commandLineArguments.useLimiter >= 0,
            "nObs"                  : commandLineArguments.useLimiter,
            "ghostLayerWidth"       : commandLineArguments.ghostLayerWidth,
            "pathToLibxsmmGemmGenerator"  : absolutePathToLibxsmm,
            "quadratureType"        : ("Gauss-Lobatto" if commandLineArguments.useGaussLobatto else "Gauss-Legendre"),
            "useCERKGuess"          : commandLineArguments.useCERKGuess,
            "useLibxsmm"            : True,
            "runtimeDebug"          : False #for debug
        }
        self.config["useSourceOrNCP"] = self.config["useSource"] or self.config["useNCP"]
        self.validateConfig(simdWidth.keys())
        self.config["vectSize"] = simdWidth[self.config["architecture"]] #only initialize once architecture has been validated
        self.baseContext = self.generateBaseContext() # default context build from config
        self.gemmList = [] #list to store the name of all generated gemms (used for gemmsCPPModel)


    def validateConfig(self, validArchitectures):
        """Ensure the configuration fit some constraint, raise errors if not"""
        if not (self.config["architecture"] in validArchitectures):
           raise ValueError("Architecture not recognized. Available architecture: "+str(validArchitectures))
        if not (self.config["numerics"] == "linear" or self.config["numerics"] == "nonlinear"):
            raise ValueError("numerics has to be linear or nonlinear")
        if self.config["nVar"] < 0:
           raise ValueError("Number of variables must be >=0 ")
        if self.config["nPar"] < 0:
           raise ValueError("Number of parameters must be >= 0")
        if self.config["nDim"] < 2 or self.config["nDim"] > 3:
           raise ValueError("Number of dimensions must be 2 or 3")
        if self.config["nDof"] < 0 or self.config["nDof"] > 9:
           raise ValueError("Order has to be between 0 and 9")
        if (self.config["useSource"] and not self.config["useSourceVect"] and self.config["useNCPVect"]) or (self.config["useNCP"] and not self.config["useNCPVect"] and self.config["useSourceVect"]) :
            raise ValueError("If using source and NCP, both or neither must be vectorized")


    def printConfig(self):
        print(self.config)


    def generateBaseContext(self):
        """Generate a base context for the models from the config (use hard copy)"""
        context = copy.copy(self.config)
        context["nVarPad"]  = self.getSizeWithPadding(context["nVar"])
        context["nParPad"]  = self.getSizeWithPadding(context["nPar"])
        context["nDataPad"] = self.getSizeWithPadding(context["nData"])
        context["nDofPad"]  = self.getSizeWithPadding(context["nDof"])
        context["nDof3D"]   = 1 if context["nDim"] == 2 else context["nDof"]
        context["isLinear"] = context["numerics"] == "linear"
        context["solverHeader"]      = context["solverName"].split("::")[1] + ".h"
        context["codeNamespaceList"] = context["codeNamespace"].split("::")
        context["guardNamespace"]    = "_".join(context["codeNamespaceList"]).upper()
        context["nDofLim"] = 2*context["nDof"]-1 #for limiter
        context["nDofLimPad"] = self.getSizeWithPadding(context["nDofLim"])
        context["nDofLim3D"] = 1 if context["nDim"] == 2 else context["nDofLim"]
        context["ghostLayerWidth3D"] = 0 if context["nDim"] == 2 else context["ghostLayerWidth"]
        context["useVectPDEs"] = context["useFluxVect"] or True #TODO JMG add other vect
        return context

    def getSizeWithPadding(self, sizeWithoutPadding):
        """Return the size of the input with the architecture specific padding added"""
        return self.config["vectSize"] * int((sizeWithoutPadding+(self.config["vectSize"]-1))/self.config["vectSize"])


    def getPadSize(self, sizeWithoutPadding):
        """Return the size of padding required for its input"""
        return self.getSizeWithPadding(sizeWithoutPadding) - sizeWithoutPadding


    def generateCode(self):
        """Main method: call the models to generate the code"""
        
        # create directory for output files if not existing
        try:
            os.makedirs(self.config['pathToOutputDirectory'])
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        
        # remove all .cpp, .cpph, .c and .h files (we are in append mode!)
        for fileName in os.listdir(self.config['pathToOutputDirectory']):
            _ , ext = os.path.splitext(fileName)
            if(ext in [".cpp", ".cpph", ".c", ".h"]):
                os.remove(self.config['pathToOutputDirectory'] + "/" + fileName)
        
        # generate new files
        runtimes = {}
        
        start = time.perf_counter()
        adjustSolution = adjustSolutionModel.AdjustSolutionModel(self.baseContext)
        adjustSolution.generateCode()
        runtimes["adjustSolution"] = time.perf_counter() - start
        
        start = time.perf_counter()
        amrRoutines = amrRoutinesModel.AMRRoutinesModel(self.baseContext, self)
        amrRoutines.generateCode()
        runtimes["amrRoutines"] = time.perf_counter() - start
        
        start = time.perf_counter()
        boundaryConditions = boundaryConditionsModel.BoundaryConditionsModel(self.baseContext)
        boundaryConditions.generateCode()
        runtimes["boundaryConditions"] = time.perf_counter() - start
        
        start = time.perf_counter()
        configurationParameters = configurationParametersModel.ConfigurationParametersModel(self.baseContext)
        configurationParameters.generateCode()
        runtimes["configurationParameters"] = time.perf_counter() - start
        
        start = time.perf_counter()
        converter = converterModel.ConverterModel(self.baseContext)
        converter.generateCode()
        runtimes["converter"] = time.perf_counter() - start
        
        start = time.perf_counter()
        deltaDistribution = deltaDistributionModel.DeltaDistributionModel(self.baseContext)
        deltaDistribution.generateCode()
        runtimes["deltaDistribution"] = time.perf_counter() - start
        
        start = time.perf_counter()
        dgMatrix = dgMatrixModel.DGMatrixModel(self.baseContext)
        dgMatrix.generateCode()
        runtimes["dgMatrix"] = time.perf_counter() - start
        
        start = time.perf_counter()
        faceIntegral = faceIntegralModel.FaceIntegralModel(self.baseContext)
        faceIntegral.generateCode()
        runtimes["faceIntegral"] = time.perf_counter() - start
        
        start = time.perf_counter()
        fusedSpaceTimePredictorVolumeIntegral = fusedSpaceTimePredictorVolumeIntegralModel.FusedSpaceTimePredictorVolumeIntegralModel(self.baseContext, self)
        fusedSpaceTimePredictorVolumeIntegral.generateCode()
        runtimes["fusedSpaceTimePredictorVolumeIntegral"] = time.perf_counter() - start
        
        start = time.perf_counter()
        kernelsHeader = kernelsHeaderModel.KernelsHeaderModel(self.baseContext)
        kernelsHeader.generateCode()
        runtimes["kernelsHeader"] = time.perf_counter() - start
        
        start = time.perf_counter()
        limiter = limiterModel.LimiterModel(self.baseContext, self)
        limiter.generateCode()
        runtimes["limiter"] = time.perf_counter() - start
        
        start = time.perf_counter()
        matrixUtils = matrixUtilsModel.MatrixUtilsModel(self.baseContext)
        matrixUtils.generateCode()
        runtimes["matrixUtils"] = time.perf_counter() - start
        
        start = time.perf_counter()
        quadrature = quadratureModel.QuadratureModel(self.baseContext, self)
        quadrature.generateCode()
        runtimes["quadrature"] = time.perf_counter() - start
        
        start = time.perf_counter()
        riemann = riemannModel.RiemannModel(self.baseContext)
        riemann.generateCode()
        runtimes["riemann"] = time.perf_counter() - start
        
        start = time.perf_counter()
        solutionUpdate = solutionUpdateModel.SolutionUpdateModel(self.baseContext)
        solutionUpdate.generateCode()
        runtimes["solutionUpdate"] = time.perf_counter() - start
        
        start = time.perf_counter()
        stableTimeStepSize = stableTimeStepSizeModel.StableTimeStepSizeModel(self.baseContext)
        stableTimeStepSize.generateCode()
        runtimes["stableTimeStepSize"] = time.perf_counter() - start
        
        start = time.perf_counter()
        surfaceIntegral = surfaceIntegralModel.SurfaceIntegralModel(self.baseContext)
        surfaceIntegral.generateCode()
        runtimes["surfaceIntegral"] = time.perf_counter() - start
        
        # must be run only after all gemm have been generated
        start = time.perf_counter()
        gemmsContext = copy.copy(self.baseContext)
        gemmsContext["gemmList"] = self.gemmList
        gemmsCPP = gemmsCPPModel.GemmsCPPModel(gemmsContext)
        gemmsCPP.generateCode()
        runtimes["gemmsCPP"] = time.perf_counter() - start
        
        if self.config['runtimeDebug']:
            for key, value in runtimes.items():
                print(key+": "+str(value))


    def generateGemms(self, outputFileName, matmulConfigList):
        """Generate the gemms with the given config using LIBXSMM"""
        for matmul in matmulConfigList:
            # add the gemm name to the list of generated gemm
            self.gemmList.append(matmul.baseroutinename)
            # for plain assembly code (rather than inline assembly) choose dense_asm
            commandLineArguments = " " + "dense"  + \
                " " + os.path.join(self.config["pathToOutputDirectory"], outputFileName) + \
                " " + self.config["codeNamespace"] + "::" + matmul.baseroutinename + \
                " " + str(matmul.M) + \
                " " + str(matmul.N) + \
                " " + str(matmul.K) + \
                " " + str(matmul.LDA) + \
                " " + str(matmul.LDB) + \
                " " + str(matmul.LDC) + \
                " " + str(matmul.alpha) + \
                " " + str(matmul.beta) + \
                " " + str(matmul.alignment_A) + \
                " " + str(matmul.alignment_C) + \
                " " + self.config["architecture"] + \
                " " + matmul.prefetchStrategy + \
                " " + "DP" #always use double precision, "SP" for single
            bashCommand = self.config["pathToLibxsmmGemmGenerator"] + commandLineArguments
            subprocess.call(bashCommand.split())
