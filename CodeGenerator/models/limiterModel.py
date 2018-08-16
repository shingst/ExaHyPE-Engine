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
# Generate the converter, used to mix generic and optimized kernels
#


from .abstractModelBaseClass import AbstractModelBaseClass

from utils.MatmulConfig import MatmulConfig


class LimiterModel(AbstractModelBaseClass):
    
    def generateCode(self):
        if(not self.context['useLimiter']):
            return None
        self.context["gemm_dg2fv"]  = "gemm_"+str(self.context["nVar"])+"_"+str(self.context["nDofLim"])+"_"+str(self.context["nDof"])   +"_dg2fv"
        self.context["gemm_fv2dg"]  = "gemm_"+str(self.context["nVar"])+"_"+str(self.context["nDof"])   +"_"+str(self.context["nDofLim"])+"_fv2dg"
        self.context["gemm_uh2lob"] = "gemm_"+str(self.context["nVar"])+"_"+str(self.context["nDof"])   +"_"+str(self.context["nDof"])   +"_uh2lob"
        
        self.render("limiter_cpp.template", "limiter.cpp")
        # generates gemms
        if(self.context["useLibxsmm"]):
            self.controller.generateGemms("asm_limiter.c", self.context["gemmList"].values())
    
    
    def buildGemmsConfig(self):
        # define a sequence of matmul configs
        self.context["gemmList"] = {}

        #-----------------------------
        # implementation file
        #-----------------------------        
        self.context["gemmList"]["dg2fv"] = MatmulConfig(  
                                    # M
                                    self.context["nVar"],      \
                                    # N
                                    self.context["nDofLim"],       \
                                    # K
                                    self.context["nDof"],       \
                                    # LDA
                                    self.context["nData"],   \
                                    # LDB
                                    self.context["nDofPad"],    \
                                    # LDC
                                    self.context["nVar"],   \
                                    # alpha 
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    1,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "dg2fv",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
        self.context["gemmList"]["fv2dg"] = MatmulConfig(  
                                    # M
                                    self.context["nVar"],      \
                                    # N
                                    self.context["nDof"],       \
                                    # K
                                    self.context["nDofLim"],       \
                                    # LDA
                                    self.context["nVar"],   \
                                    # LDB
                                    self.context["nDofLimPad"],    \
                                    # LDC
                                    self.context["nData"],   \
                                    # alpha 
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    1,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "fv2dg",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
        self.context["gemmList"]["uh2lob"] = MatmulConfig(  
                                    # M
                                    self.context["nVar"],      \
                                    # N
                                    self.context["nDof"],       \
                                    # K
                                    self.context["nDof"],       \
                                    # LDA
                                    self.context["nData"],   \
                                    # LDB
                                    self.context["nDofPad"],    \
                                    # LDC
                                    self.context["nVar"],   \
                                    # alpha 
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    1,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "uh2lob",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
