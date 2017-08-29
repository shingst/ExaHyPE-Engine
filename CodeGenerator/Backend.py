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
# This file is pivotal to the code generator. It manages 
# the internal decisions about padding, triggers the code
# generation of the solver kernels and calls the back end 
# assembly code generation.
#

import os
from os.path import join
from os.path import isfile
import copy
import subprocess
import errno
from glob import iglob
from shutil import move
import KernelsHeaderGenerator
import SpaceTimePredictorGenerator
import VolumeIntegralGenerator
import SurfaceIntegralGenerator
import RiemannGenerator
import SolutionUpdateGenerator
import AdjustSolutionGenerator
import StableTimeStepSizeGenerator
import WeightsGenerator
import DGMatrixGenerator
import ConfigurationParametersGenerator
import BoundaryConditionsGenerator
import ConverterGenerator
import string
import re

m_architecture           = ''
m_precision              = ''
m_config                 = {}
m_numerics               = ''
m_pathToLibxsmmGenerator = ''
m_simdWidth              =  {'SP':  {'noarch' : 1,
                                     'wsm'    : 4,
                                     'snb'    : 8,
                                     'hsw'    : 8,
                                     'knc'    : 16,
                                     'knl'    : 16 },
                             'DP': { 'noarch' : 1,
                                     'wsm'    : 2,
                                     'snb'    : 4,
                                     'hsw'    : 4,
                                     'knc'    : 8,
                                     'knl'    : 8 }
                            }


def executeLibxsmmGenerator(i_commandLineParameters):
    l_bashCommand = m_pathToLibxsmmGenerator + '/libxsmm_gemm_generator ' + i_commandLineParameters
    subprocess.call(l_bashCommand.split())


def generateAssemblerCode(i_pathToOutputFile,
                          i_matmulConfigList):
    l_pathToAsmFile = os.path.splitext(i_pathToOutputFile)[0]+'.c'
    for l_matmul in i_matmulConfigList:
        # for plain assembly code (rather than inline assembly) choose dense_asm
        l_commandLineArguments =       "dense"  + \
                                 ' ' + m_pathToLibxsmmGenerator+"/"+l_pathToAsmFile + \
                                 ' ' + m_config['codeNamespace'] + '::' + l_matmul.baseroutinename + \
                                 ' ' + str(l_matmul.M) + \
                                 ' ' + str(l_matmul.N) + \
                                 ' ' + str(l_matmul.K) + \
                                 ' ' + str(l_matmul.LDA) + \
                                 ' ' + str(l_matmul.LDB) + \
                                 ' ' + str(l_matmul.LDC) + \
                                 ' ' + str(l_matmul.alpha) + \
                                 ' ' + str(l_matmul.beta) + \
                                 ' ' + str(l_matmul.alignment_A) + \
                                 ' ' + str(l_matmul.alignment_C) + \
                                 ' ' + m_architecture + \
                                 ' ' + l_matmul.prefetchStrategy+ \
                                 ' ' + m_precision 
        executeLibxsmmGenerator(l_commandLineArguments)


def getSizeWithPadding(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = m_simdWidth[m_precision][m_architecture]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    return l_sizeWithPadding


def getPadWidth(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = m_simdWidth[m_precision][m_architecture]
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
        if(l_ext in ['.cpp', '.cpph', '.c', '.h']):
            os.remove(i_outputDirectory + "/" + l_fileName)


def executeBashCommand(i_command, i_commandLineParameters):
    # usage: executeBashCommand("ls", "-l -a")
    l_bashCommand = i_command + " " + i_commandLineParameters
    l_commandOutput = subprocess.check_output(l_bashCommand, shell=True)
    return l_commandOutput


def generateContext(i_config):
    context = copy.copy(i_config)
    context['nVarPad'] = getSizeWithPadding(context['nVar'])
    context['nDofPad'] = getSizeWithPadding(context['nDof'])
    context['nDof3D'] = 1 if context['nDim'] == 2 else context['nDof']
    context['isLinear'] = context['numerics'] == "linear"
    context['solverHeader'] = context['solverName'].split('::')[1] + '.h'
    #context['FloatingPointFormat'] = 'float' if 'm_precision' == 'SP' else 'double'
    context['codeNamespaceList'] = context['codeNamespace'].split('::')
    context['guardNamespace'] = '_'.join(context['codeNamespaceList']).upper()
    return context

    
def generateComputeKernels():
    kernelsHeaderGenerator = KernelsHeaderGenerator.KernelsHeaderGenerator(generateContext(m_config))
    kernelsHeaderGenerator.generateCode()
    spaceTimePredictorGenerator = SpaceTimePredictorGenerator.SpaceTimePredictorGenerator(generateContext(m_config))
    spaceTimePredictorGenerator.generateCode()
    volumeIntegralGenerator = VolumeIntegralGenerator.VolumeIntegralGenerator(generateContext(m_config))
    volumeIntegralGenerator.generateCode()
    surfaceIntegralGenerator = SurfaceIntegralGenerator.SurfaceIntegralGenerator(generateContext(m_config))
    surfaceIntegralGenerator.generateCode()
    riemannGenerator = RiemannGenerator.RiemannGenerator(generateContext(m_config))
    riemannGenerator.generateCode()
    solutionUpdateGenerator = SolutionUpdateGenerator.SolutionUpdateGenerator(generateContext(m_config))
    solutionUpdateGenerator.generateCode()
    adjustSolutionGenerator = AdjustSolutionGenerator.AdjustSolutionGenerator(generateContext(m_config))
    adjustSolutionGenerator.generateCode()
    stableTimeStepSizeGenerator = StableTimeStepSizeGenerator.StableTimeStepSizeGenerator(generateContext(m_config))
    stableTimeStepSizeGenerator.generateCode()
    weightsGenerator = WeightsGenerator.WeightsGenerator(generateContext(m_config))
    weightsGenerator.generateCode()
    dgMatrixGenerator = DGMatrixGenerator.DGMatrixGenerator(generateContext(m_config))
    dgMatrixGenerator.generateCode()
    # no ccph anymore => not needed anymore. Legacy TODO JMG clean later
    #cpphGemmsGenerator = CpphGemmsGenerator.CpphGemmsGenerator(generateContext(m_config))
    #cpphGemmsGenerator.generateCode()
    configurationParametersGenerator = ConfigurationParametersGenerator.ConfigurationParametersGenerator(generateContext(m_config))
    configurationParametersGenerator.generateCode()
    boundaryConditionsGenerator = BoundaryConditionsGenerator.BoundaryConditionsGenerator(generateContext(m_config))
    boundaryConditionsGenerator.generateCode()
    converterGenerator = ConverterGenerator.ConverterGenerator(generateContext(m_config))
    converterGenerator.generateCode()



def moveGeneratedFiles(i_pathToSrc,i_pathToDest):
    l_fileTypes = ('*.h', '*.cpp', '*.c', '*.cpph')
    l_fileList = []
    for l_file in l_fileTypes:
        l_fileList.extend(iglob(i_pathToSrc+"/"+l_file))

    for l_file in l_fileList:
        if(isfile(l_file)):
            move(l_file, i_pathToDest)

# -------------------------------------------------------------------
# Class variables 
# -------------------------------------------------------------------
def setArchitecture(i_architecture):
    global m_architecture 
    m_architecture = i_architecture

def setPrecision(i_precision):
    global m_precision
    m_precision = i_precision

def setConfig(i_config):
    global m_config
    m_config = i_config

def setNumerics(i_numerics):
    global m_numerics
    m_numerics = i_numerics

def setPathToLibxsmmGenerator(i_pathToLibxsmmGenerator):
    global m_pathToLibxsmmGenerator
    m_pathToLibxsmmGenerator = i_pathToLibxsmmGenerator

# -------------------------------------------------------------------
# helpers
# -------------------------------------------------------------------
def reindentLine(i_line, i_nSpaces=8):
    return (' ' * i_nSpaces) + i_line

def reindentBlock(i_string, i_nSpaces):
    l_stringList = i_string.split('\n')
    l_stringList = list(map(lambda line, ns=i_nSpaces: reindentLine(line, ns), l_stringList))
    l_string = "\n".join(l_stringList)
    return l_string
