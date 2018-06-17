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
# Generates the code for the AMR
#


from .abstractModelBaseClass import AbstractModelBaseClass

import controller
from utils.MatmulConfig import MatmulConfig


class AMRRoutinesModel(AbstractModelBaseClass):    
    
    def generateCode(self):
        self.context["gemm_face_Q"] = "gemm_"+str(self.context["nDataPad"])+"_"+str(self.context["nDof"])+"_"+str(self.context["nDof"])+"_face_Q"
        self.context["gemm_face_F"] = "gemm_"+str(self.context["nVarPad"]) +"_"+str(self.context["nDof"])+"_"+str(self.context["nDof"])+"_face_F"
        self.context["gemm_volume"] = "gemm_"+str(self.context["nData"])+"_"+str(self.context["nDof"])+"_"+str(self.context["nDof"])+"_volume"
        
        self.render("amrRoutines_cpp.template", "amrRoutines.cpp")
        # generates gemms
        if(self.context["useLibxsmm"]):
            self.controller.generateGemms("asm_amrRoutines.c", self.context["gemmList"].values())
    
    
    def buildGemmsConfig(self):
        # define a sequence of matmul configs
        self.context["gemmList"] = {}

        #-----------------------------
        # implementation file
        #-----------------------------
        self.context["gemmList"]["face_Q"] = MatmulConfig(  
                                    # M
                                    self.context["nDataPad"],      \
                                    # N
                                    self.context["nDof"],       \
                                    # K
                                    self.context["nDof"],       \
                                    # LDA
                                    self.context["nDataPad"],   \
                                    # LDB
                                    self.context["nDofPad"],    \
                                    # LDC
                                    self.context["nDataPad"],   \
                                    # alpha
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    0,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "face_Q",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
        self.context["gemmList"]["face_F"] = MatmulConfig(  
                                    # M
                                    self.context["nVarPad"],       \
                                    # N
                                    self.context["nDof"],       \
                                    # K
                                    self.context["nDof"],       \
                                    # LDA
                                    self.context["nVarPad"],    \
                                    # LDB
                                    self.context["nDofPad"],    \
                                    # LDC
                                    self.context["nVarPad"],    \
                                    # alpha
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    0,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "face_F",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
        self.context["gemmList"]["volume"] = MatmulConfig(  
                                    # M
                                    self.context["nData"],       \
                                    # N
                                    self.context["nDof"],       \
                                    # K
                                    self.context["nDof"],       \
                                    # LDA
                                    self.context["nData"],      \
                                    # LDB
                                    self.context["nDofPad"],    \
                                    # LDC
                                    self.context["nData"],      \
                                    # alpha
                                    1,                            \
                                    # beta
                                    1,                            \
                                    # alignment A
                                    0,                            \
                                    # alignment C
                                    1,                            \
                                    # name
                                    "volume",                     \
                                    # prefetching
                                    "nopf",                       \
                                    # type
                                    "gemm")
