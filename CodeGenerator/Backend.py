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
# This file is pivotal to the code generator. It manages 
# the internal decisions about padding, triggers the code
# generation of the solver kernels and calls the back end 
# assembly code generation.
#


import os
import copy
import subprocess
import errno
import time

from generators import *


g_runtimes               = {}
g_config                 = {}
g_simdWidth              = {  "noarch" : 1,
                              "wsm"    : 2,
                              "snb"    : 4,
                              "hsw"    : 4,
                              "knc"    : 8,
                              "knl"    : 8 
                            }

def setConfig(i_config):
    global g_config
    g_config = i_config

    
def getSizeWithPadding(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = g_simdWidth[g_config["architecture"]]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    return l_sizeWithPadding


def getPadWidth(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = g_simdWidth[g_config["architecture"]]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    l_padWidth        = l_sizeWithPadding - i_sizeWithoutPadding
    return l_padWidth


def prepareOutputDirectory(i_outputDirectory):
    # create directory for output files if not existing
    try:
        os.makedirs(i_outputDirectory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    # remove all .cpp, .cpph, .c and .h files (we are in append mode!)
    for l_fileName in os.listdir(i_outputDirectory):
        _ , l_ext = os.path.splitext(l_fileName)
        if(l_ext in [".cpp", ".cpph", ".c", ".h"]):
            os.remove(i_outputDirectory + "/" + l_fileName)


def generateContext(i_config):
    context = copy.copy(i_config)
    context["nVarPad"]  = getSizeWithPadding(context["nVar"])
    context["nParPad"]  = getSizeWithPadding(context["nPar"])
    context["nDataPad"] = getSizeWithPadding(context["nData"])
    context["nDofPad"]  = getSizeWithPadding(context["nDof"])
    context["nDof3D"]   = 1 if context["nDim"] == 2 else context["nDof"]
    context["isLinear"] = context["numerics"] == "linear"
    context["solverHeader"]      = context["solverName"].split("::")[1] + ".h"
    context["codeNamespaceList"] = context["codeNamespace"].split("::")
    context["guardNamespace"]    = "_".join(context["codeNamespaceList"]).upper()
    context["nDofLim"] = 2*context["nDof"]-1 #for limiter
    context["nDofLimPad"] = getSizeWithPadding(context["nDofLim"])
    context["nDofLim3D"] = 1 if context["nDim"] == 2 else context["nDofLim"]
    context["ghostLayerWidth3D"] = 0 if context["nDim"] == 2 else context["ghostLayerWidth"]
    return context

    
def generateComputeKernels():
    start = time.perf_counter()
    kernelsHeaderGenerator = KernelsHeaderGenerator.KernelsHeaderGenerator(generateContext(g_config))
    kernelsHeaderGenerator.generateCode()
    g_runtimes["kernelsHeaderGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    fusedSpaceTimePredictorVolumeIntegralGenerator = FusedSpaceTimePredictorVolumeIntegralGenerator.FusedSpaceTimePredictorVolumeIntegralGenerator(generateContext(g_config))
    fusedSpaceTimePredictorVolumeIntegralGenerator.generateCode()
    g_runtimes["fusedSpaceTimePredictorVolumeIntegralGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    surfaceIntegralGenerator = SurfaceIntegralGenerator.SurfaceIntegralGenerator(generateContext(g_config))
    surfaceIntegralGenerator.generateCode()
    g_runtimes["surfaceIntegralGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    riemannGenerator = RiemannGenerator.RiemannGenerator(generateContext(g_config))
    riemannGenerator.generateCode()
    g_runtimes["riemannGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    solutionUpdateGenerator = SolutionUpdateGenerator.SolutionUpdateGenerator(generateContext(g_config))
    solutionUpdateGenerator.generateCode()
    g_runtimes["solutionUpdateGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    adjustSolutionGenerator = AdjustSolutionGenerator.AdjustSolutionGenerator(generateContext(g_config))
    adjustSolutionGenerator.generateCode()
    g_runtimes["adjustSolutionGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    stableTimeStepSizeGenerator = StableTimeStepSizeGenerator.StableTimeStepSizeGenerator(generateContext(g_config))
    stableTimeStepSizeGenerator.generateCode()
    g_runtimes["stableTimeStepSizeGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    quadratureGenerator = QuadratureGenerator.QuadratureGenerator(generateContext(g_config))
    quadratureGenerator.generateCode()
    g_runtimes["quadratureGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    dgMatrixGenerator = DGMatrixGenerator.DGMatrixGenerator(generateContext(g_config))
    dgMatrixGenerator.generateCode()
    g_runtimes["dgMatrixGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    gemmsCPPGenerator = GemmsCPPGenerator.GemmsCPPGenerator(generateContext(g_config))
    gemmsCPPGenerator.generateCode()
    g_runtimes["gemmsCPPGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    configurationParametersGenerator = ConfigurationParametersGenerator.ConfigurationParametersGenerator(generateContext(g_config))
    configurationParametersGenerator.generateCode()
    g_runtimes["configurationParametersGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    boundaryConditionsGenerator = BoundaryConditionsGenerator.BoundaryConditionsGenerator(generateContext(g_config))
    boundaryConditionsGenerator.generateCode()
    g_runtimes["boundaryConditionsGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    converterGenerator = ConverterGenerator.ConverterGenerator(generateContext(g_config))
    converterGenerator.generateCode()
    g_runtimes["converterGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    amrRoutinesGenerator = AMRRoutinesGenerator.AMRRoutinesGenerator(generateContext(g_config))
    amrRoutinesGenerator.generateCode()
    g_runtimes["amrRoutinesGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    deltaDistributionGenerator = DeltaDistributionGenerator.DeltaDistributionGenerator(generateContext(g_config))
    deltaDistributionGenerator.generateCode()
    g_runtimes["deltaDistributionGenerator"] = time.perf_counter() - start
    start = time.perf_counter()
    limiterGenerator = LimiterGenerator.LimiterGenerator(generateContext(g_config))
    limiterGenerator.generateCode()
    g_runtimes["LimiterGenerator"] = time.perf_counter() - start


def executeLibxsmmGenerator(i_commandLineParameters):
    l_bashCommand = g_config["pathToLibxsmmGemmGenerator"] + i_commandLineParameters
    subprocess.call(l_bashCommand.split())


def generateAssemblerCode(i_outputFileName,
                          i_matmulConfigList):
    l_pathToAsmFile = os.path.splitext(i_outputFileName)[0]+".c"
    for l_matmul in i_matmulConfigList:
        # for plain assembly code (rather than inline assembly) choose dense_asm
        l_commandLineArguments = " " + "dense"  + \
                                 " " + os.path.join(g_config["pathToOutputDirectory"],l_pathToAsmFile) + \
                                 " " + g_config["codeNamespace"] + "::" + l_matmul.baseroutinename + \
                                 " " + str(l_matmul.M) + \
                                 " " + str(l_matmul.N) + \
                                 " " + str(l_matmul.K) + \
                                 " " + str(l_matmul.LDA) + \
                                 " " + str(l_matmul.LDB) + \
                                 " " + str(l_matmul.LDC) + \
                                 " " + str(l_matmul.alpha) + \
                                 " " + str(l_matmul.beta) + \
                                 " " + str(l_matmul.alignment_A) + \
                                 " " + str(l_matmul.alignment_C) + \
                                 " " + g_config["architecture"] + \
                                 " " + l_matmul.prefetchStrategy+ \
                                 " " + "DP" #always use double precision, "SP" for single
        executeLibxsmmGenerator(l_commandLineArguments)

def printRuntimes():
    for key, value in g_runtimes.items():
        print(key+": "+str(value))
