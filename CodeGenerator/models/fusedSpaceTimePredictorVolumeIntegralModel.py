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
# Generates the SpaceTimePredictor Kernel
#
# Call user function flux, source, ncp
#


from .abstractModelBaseClass import AbstractModelBaseClass

import copy

import controller
from utils.MatmulConfig import MatmulConfig


class FusedSpaceTimePredictorVolumeIntegralModel(AbstractModelBaseClass):  
    
    def generateCode(self):
        gemmName = "gemm_"+str(self.context["nVar"])+"_"+str(self.context["nDof"])+"_"+str(self.context["nDof"])
        gemmNamePad = "gemm_"+str(self.context["nVarPad"])+"_"+str(self.context["nDof"])+"_"+str(self.context["nDof"])
        self.context["gemm_gradQ_x"] = gemmName+"_gradQ_x"
        self.context["gemm_gradQ_y"] = gemmName+"_gradQ_y"
        self.context["gemm_gradQ_z"] = gemmName+"_gradQ_z"
        self.context["nVarMinusOne_seq"] = range(self.context["nVar"] - 1)
        
        if(self.context["isLinear"]):
            self.context["ncpOutputShift"] = self.controller.getSizeWithPadding(self.context["nVar"]*self.context["nDim"]) #shift used to split the tmpArray into input and output for NCP
            # size of the tmpArray
            self.context["tmpArraySize"] = max((self.context["nDof"]*self.context["nVarPad"] if self.context["useFlux"]          else 0), \
                                                 (self.context["nVar"]*self.context["nDim"]    if self.context["useMaterialParam"] else 0), \
                                                 (2*self.context["ncpOutputShift"]               if self.context["useNCP"]           else 0))
            self.context["gemm_flux_x"] = gemmNamePad+"_flux_x"
            self.context["gemm_flux_y"] = gemmNamePad+"_flux_y"
            self.context["gemm_flux_z"] = gemmNamePad+"_flux_z"
     
            self.render("fusedSPTVI_linear_cpp.template", "fusedSpaceTimePredictorVolumeIntegral.cpp")
            
            if(self.context["usePointSources"]):
                localContext = copy.copy(self.context)
                localContext["usePointSources"] = False
                localContext["nameSuffix"] = "_WithoutPS"
                
                self.render("fusedSPTVI_linear_cpp.template", "fusedSpaceTimePredictorVolumeIntegral_WithoutPS.cpp", localContext)
                
        else:
            self.context["nDof_seq"] = range(0,self.context["nDof"])
            self.context["gemm_rhs_x"] = gemmNamePad+"_rhs_x"
            self.context["gemm_rhs_y"] = gemmNamePad+"_rhs_y"
            self.context["gemm_rhs_z"] = gemmNamePad+"_rhs_z"
            self.context["gemm_gradF_x"] = gemmName+"_gradF_x"
            self.context["gemm_gradF_y"] = gemmName+"_gradF_y"
            self.context["gemm_gradF_z"] = gemmName+"_gradF_z"
            self.context["gemm_lqi"]   = gemmName+"_lqi"
            self.context["gemm_x"] = gemmNamePad+"_lduh_x"
            self.context["gemm_y"] = gemmNamePad+"_lduh_y"
            self.context["gemm_z"] = gemmNamePad+"_lduh_z"             
            self.context["i_seq"] = range(0,self.context["nDof"])
            self.context["j_seq"] = range(0,self.context["nDof"]) if (self.context["nDim"] >= 3) else [0]
            
            self.render("fusedSPTVI_nonlinear_cpp.template", "fusedSpaceTimePredictorVolumeIntegral.cpp")
       
        # generates gemms
        if(self.context["useLibxsmm"]):
            self.controller.generateGemms("asm_fstpvi.c", self.context["gemmList"].values())
    
    
    def buildGemmsConfig(self):
        self.context["gemmList"] = {}
        
        if(self.context["isLinear"]):
            if(self.context["useFlux"]):
                self.context["gemmList"]["flux_x"] = MatmulConfig(    
                                            # M
                                            self.context["nVarPad"],   \
                                            # N
                                            self.context["nDof"],      \
                                            # K
                                            self.context["nDof"],      \
                                            # LDA
                                            self.context["nVarPad"],   \
                                            # LDB
                                            self.context["nDofPad"],   \
                                            # LDC
                                            self.context["nVarPad"],   \
                                            # alpha
                                            1,                         \
                                            # beta, 0 => overwrite C
                                            0,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "flux_x",                  \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                self.context["gemmList"]["flux_y"] = MatmulConfig(    
                                            # M
                                            self.context["nVarPad"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nVarPad"] * self.context["nDof"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"], \
                                            # alpha
                                            1,                         \
                                            # beta, 0 => overwrite C
                                            0,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "flux_y",                  \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["flux_z"] = MatmulConfig(    
                                                # M
                                                self.context["nVarPad"],    \
                                                # N
                                                self.context["nDof"],    \
                                                # K
                                                self.context["nDof"],    \
                                                # LDA
                                                self.context["nVarPad"] * (self.context["nDof"]**2), \
                                                # LDB
                                                self.context["nDofPad"], \
                                                # LDC
                                                self.context["nVarPad"], \
                                                # alpha
                                                1,                         \
                                                # beta, 0 => overwrite C
                                                0,                         \
                                                # alignment A
                                                1,                         \
                                                # alignment C
                                                1,                         \
                                                # name
                                                "flux_z",                  \
                                                # prefetching
                                                "nopf",                    \
                                                # type
                                                "gemm")
            if(self.context["useNCP"]):
                self.context["gemmList"]["gradQ_x"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nDataPad"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradQ_x",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                self.context["gemmList"]["gradQ_y"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nDataPad"] * self.context["nDof"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"] * self.context["nDof"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradQ_y",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["gradQ_z"] = MatmulConfig(    
                                                # M
                                                self.context["nVar"],    \
                                                # N
                                                self.context["nDof"],    \
                                                # K
                                                self.context["nDof"],    \
                                                # LDA
                                                self.context["nDataPad"] * (self.context["nDof"] ** 2), \
                                                # LDB
                                                self.context["nDofPad"], \
                                                # LDC
                                                self.context["nVarPad"] * (self.context["nDof"] ** 2), \
                                                # alpha
                                                1,                         \
                                                # beta
                                                1,                         \
                                                # alignment A
                                                1,                         \
                                                # alignment C
                                                1,                         \
                                                # name
                                                "gradQ_z",                   \
                                                # prefetching
                                                "nopf",                    \
                                                # type
                                                "gemm")
        else: #NonLinear
            if(self.context["useFlux"]):
                self.context["gemmList"]["rhs_x"] = MatmulConfig(    
                                            # M
                                            self.context["nVarPad"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nVarPad"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "rhs_x",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                self.context["gemmList"]["rhs_y"] = MatmulConfig(    
                                            # M
                                            self.context["nVarPad"],                             \
                                            # N
                                            self.context["nDof"],                             \
                                            # K
                                            self.context["nDof"],                             \
                                            # LDA
                                            self.context["nVarPad"]* self.context["nDof"],     \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"] * self.context["nDof"],     \
                                            # alpha
                                            1,                                                 \
                                            # beta
                                            1,                                                 \
                                            # alignment A
                                            1,                                                 \
                                            # alignment C
                                            1,                                                 \
                                            # name
                                            "rhs_y",                                           \
                                            # prefetching
                                            "nopf",                                            \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["rhs_z"] = MatmulConfig(    
                                                # M
                                                self.context["nVarPad"],                             \
                                                # N
                                                self.context["nDof"],                             \
                                                # K
                                                self.context["nDof"],                             \
                                                # LDA
                                                self.context["nVarPad"] * (self.context["nDof"]**2),     \
                                                # LDB
                                                self.context["nDofPad"],                          \
                                                # LDC
                                                self.context["nVarPad"] * (self.context["nDof"]**2),  \
                                                # alpha
                                                1,                                                 \
                                                # beta
                                                1,                                                 \
                                                # alignment A
                                                1,                                                 \
                                                # alignment C
                                                1,                                                 \
                                                # name
                                                "rhs_z",                                           \
                                                # prefetching
                                                "nopf",                                            \
                                                # type
                                                "gemm")
                self.context["gemmList"]["lduh_x"] = MatmulConfig(  
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
                                            self.context["nVarPad"],       \
                                            # alpha 
                                            1,                            \
                                            # beta
                                            1,                            \
                                            # alignment A
                                            0,                            \
                                            # alignment C
                                            0,                            \
                                            # name
                                            "lduh_x",                     \
                                            # prefetching
                                            "nopf",                       \
                                            # type
                                            "gemm")
                self.context["gemmList"]["lduh_y"] = MatmulConfig(  
                                            # M
                                            self.context["nVarPad"],                         \
                                            # N
                                            self.context["nDof"],                         \
                                            # K
                                            self.context["nDof"],                         \
                                            # LDA
                                            self.context["nVarPad"],                      \
                                            # LDB
                                            self.context["nDofPad"],                      \
                                            # LDC
                                            self.context["nVarPad"]*self.context["nDof"],  \
                                            # alpha 
                                            1,                                              \
                                            # beta
                                            1,                                              \
                                            # alignment A
                                            0,                                              \
                                            # alignment C
                                            0,                                              \
                                            # name
                                            "lduh_y",                                       \
                                            # prefetching
                                            "nopf",                                         \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["lduh_z"] = MatmulConfig(  
                                                # M
                                                self.context["nVarPad"],                             \
                                                # N
                                                self.context["nDof"],                             \
                                                # K
                                                self.context["nDof"],                             \
                                                # LDA
                                                self.context["nVarPad"],                          \
                                                # LDB
                                                self.context["nDofPad"],                          \
                                                # LDC
                                                self.context["nVarPad"]*(self.context["nDof"]**2), \
                                                # alpha 
                                                1,                                                  \
                                                # beta
                                                1,                                                  \
                                                # alignment A
                                                0,                                                  \
                                                # alignment C
                                                0,                                                  \
                                                # name
                                                "lduh_z",                                           \
                                                # prefetching
                                                "nopf",                                             \
                                                # type
                                                "gemm")
                self.context["gemmList"]["gradF_x"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nVarPad"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradF_x",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                self.context["gemmList"]["gradF_y"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nVarPad"] * self.context["nDof"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"] * self.context["nDof"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradF_y",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["gradF_z"] = MatmulConfig(    
                                                # M
                                                self.context["nVar"],    \
                                                # N
                                                self.context["nDof"],    \
                                                # K
                                                self.context["nDof"],    \
                                                # LDA
                                                self.context["nVarPad"] * (self.context["nDof"] ** 2), \
                                                # LDB
                                                self.context["nDofPad"], \
                                                # LDC
                                                self.context["nVarPad"] * (self.context["nDof"] ** 2), \
                                                # alpha
                                                1,                         \
                                                # beta
                                                1,                         \
                                                # alignment A
                                                1,                         \
                                                # alignment C
                                                1,                         \
                                                # name
                                                "gradF_z",                   \
                                                # prefetching
                                                "nopf",                    \
                                                # type
                                                "gemm")
            if(self.context["useNCP"]):
                self.context["gemmList"]["gradQ_x"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nDataPad"] * self.context["nDof"], \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"] * self.context["nDim"] * self.context["nDof"], \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradQ_x",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                self.context["gemmList"]["gradQ_y"] = MatmulConfig(    
                                            # M
                                            self.context["nVar"],    \
                                            # N
                                            self.context["nDof"],    \
                                            # K
                                            self.context["nDof"],    \
                                            # LDA
                                            self.context["nDataPad"] * (self.context["nDof"] ** 2), \
                                            # LDB
                                            self.context["nDofPad"], \
                                            # LDC
                                            self.context["nVarPad"] * self.context["nDim"] * (self.context["nDof"] ** 2), \
                                            # alpha
                                            1,                         \
                                            # beta
                                            1,                         \
                                            # alignment A
                                            1,                         \
                                            # alignment C
                                            1,                         \
                                            # name
                                            "gradQ_y",                   \
                                            # prefetching
                                            "nopf",                    \
                                            # type
                                            "gemm")
                if(self.context["nDim"]>=3):
                    self.context["gemmList"]["gradQ_z"] = MatmulConfig(    
                                                # M
                                                self.context["nVar"],    \
                                                # N
                                                self.context["nDof"],    \
                                                # K
                                                self.context["nDof"],    \
                                                # LDA
                                                self.context["nDataPad"] * (self.context["nDof"] ** 3), \
                                                # LDB
                                                self.context["nDofPad"], \
                                                # LDC
                                                self.context["nVarPad"] * self.context["nDim"] * (self.context["nDof"] ** 3), \
                                                # alpha
                                                1,                         \
                                                # beta
                                                1,                         \
                                                # alignment A
                                                1,                         \
                                                # alignment C
                                                1,                         \
                                                # name
                                                "gradQ_z",                   \
                                                # prefetching
                                                "nopf",                    \
                                                # type
                                                "gemm")
            self.context["gemmList"]["lqi"] = MatmulConfig(    
                                        # M
                                        self.context["nVar"],                             \
                                        # N
                                        self.context["nDof"],                             \
                                        # K
                                        self.context["nDof"],                             \
                                        # LDA
                                        self.context["nVarPad"]*(self.context["nDof"]**self.context["nDim"]), \
                                        # LDB
                                        self.context["nDofPad"], \
                                        # LDC
                                        self.context["nVarPad"], \
                                        # alpha
                                        1,                                                 \
                                        # beta
                                        0,                                                 \
                                        # alignment A
                                        1,                                                 \
                                        # alignment C
                                        1,                                                 \
                                        # name
                                        "lqi",                                             \
                                        # prefetching
                                        "nopf",                                            \
                                        # type
                                        "gemm")
