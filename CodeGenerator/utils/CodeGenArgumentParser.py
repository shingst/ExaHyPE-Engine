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
# Validates the input arguments given to the code generator.
#


import os
import AvailableConfigs


def validateLibxsmmGenerator(i_parser, i_arg):
    if not os.path.isdir(i_arg):
        i_parser.error("The libxsmm directory {} does not exist".format(i_arg))
    else:
        # directory exists, but is it a valid path to the libxsmm generator backend?
        l_pathToLibxsmm = i_arg
        l_pathToLibxsmmGenerator = l_pathToLibxsmm + "/bin/libxsmm_gemm_generator"
        if(os.path.isfile(l_pathToLibxsmmGenerator)):
            return l_pathToLibxsmm+"/bin"
        else:
            i_parser.error("Can't find the code generator of libxsmm. Did you run 'make generator' to build the library?")


def validateArchitecture(i_parser, i_arg):
    l_architecture = str(i_arg)

    # when we have to postprocess the generated assembly code we may
    # support only a subset of the available microarchitectures
    if(l_architecture not in AvailableConfigs.architectures):
        print("Driver: Unkown or unsupported microarchitecture. Continue with noarch")
        l_architecture = 'noarch'

    return l_architecture


def validatePrecision(i_parser, i_arg):
    l_precision = str(i_arg)

    if(l_precision not in AvailableConfigs.precisions):
        print("Unknown precision specified. Continue with double precision")
        l_precision    = 'DP' 

    return l_precision


def validateNumerics(i_parser, i_arg):
    l_numericsType = i_arg

    if(l_numericsType not in AvailableConfigs.numerics):
        i_parser.error("Numerics not supported. Available options are " + str(AvailableConfigs.numerics).format(i_arg))

    return l_numericsType


def validatePDE(i_parser, i_arg):
    l_pdeType = i_arg

    if(l_pdeType not in AvailableConfigs.supportedPDEs):
        i_parser.error("Driver: Unkown PDE requested. Available options are " + str(AvailableConfigs.supportedPDEs).format(i_arg))

    return l_pdeType
