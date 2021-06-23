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
# Generate the a gemms using LIBXSMM
#

import os
import subprocess

from .abstractModelBaseClass import AbstractModelBaseClass


class GemmsGeneratorModel(AbstractModelBaseClass):

    # assune context has a list of matmul config
    def generateCode(self):
        gemmList = []
        for outputFileName,matmul in self.context["matmulList"]:
            if self.context["useLibxsmm"]:
                gemmList.append(self.generateLIBXSMMgemm(outputFileName,matmul))
        return {"gemmList": gemmList}
            
    def generateLIBXSMMgemm(self, outputFileName,  matmul):
        type = "dense" # for plain assembly code (rather than inline assembly) choose dense_asm
        
        # if matmul.prefetchInput == "A":
            # prefetching = "AL2"
        # elif matmul.prefetchInput == "B":
            # prefetching = "BL2viaC"
        # elif matmul.prefetchInput == "AB":
            # prefetching = "AL2_BL2viaC"
        # else:
            # prefetching = "nopf"
        
        prefetching = "nopf"
        commandLineArguments = " " + type  + \
            " " + os.path.join(self.context["pathToOutputDirectory"], outputFileName) + \
            " " + self.context["codeNamespace"] + "::" + matmul.baseroutinename + \
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
            " " + self.context["architecture"] + \
            " " + prefetching + \
            " " + matmul.precision
        
        bashCommand = self.context["pathToLibxsmmGemmGenerator"] + commandLineArguments
        subprocess.call(bashCommand.split())
        
        return (matmul.baseroutinename, matmul.precision, prefetching != "nopf")
