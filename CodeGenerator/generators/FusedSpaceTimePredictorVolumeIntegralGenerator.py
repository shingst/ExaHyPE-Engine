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


import copy

import Backend
from utils import TemplatingUtils
from utils.MatmulConfig import MatmulConfig


class FusedSpaceTimePredictorVolumeIntegralGenerator:
    m_context = {}

    # name of generated output file
    m_filename       = "fusedSpaceTimePredictorVolumeIntegral.cpp"
    m_filename_noPS  = "fusedSpaceTimePredictorVolumeIntegral_WithoutPS.cpp"
    
    m_filename_asm   = "asm_fstpvi" 

    
    def __init__(self, i_config):
        self.m_context = i_config


    def generateCode(self):         
        gemmName = "gemm_"+str(self.m_context["nVar"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])
        gemmNamePad = "gemm_"+str(self.m_context["nVarPad"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])
        self.m_context["gemm_gradQ_x"] = gemmName+"_gradQ_x"
        self.m_context["gemm_gradQ_y"] = gemmName+"_gradQ_y"
        self.m_context["gemm_gradQ_z"] = gemmName+"_gradQ_z"
        if(self.m_context["isLinear"]):
            self.m_context["ncpOutputShift"] = Backend.getSizeWithPadding(self.m_context["nVar"]*self.m_context["nDim"]) #shift used to split the tmpArray into input and output for NCP
            # size of the tmpArray
            self.m_context["tmpArraySize"] = max((self.m_context["nDof"]*self.m_context["nVarPad"] if self.m_context["useFlux"]          else 0), \
                                                 (self.m_context["nVar"]*self.m_context["nDim"]    if self.m_context["useMaterialParam"] else 0), \
                                                 (2*self.m_context["ncpOutputShift"]               if self.m_context["useNCP"]           else 0))
            self.m_context["gemm_flux_x"] = gemmName+"_flux_x"
            self.m_context["gemm_flux_y"] = gemmName+"_flux_y"
            self.m_context["gemm_flux_z"] = gemmName+"_flux_z"
            TemplatingUtils.renderAsFile("fusedSPTVI_linear_cpp.template", self.m_filename, self.m_context)
            if(self.m_context["usePointSources"]):
                localContext = copy.copy(self.m_context)
                localContext["usePointSources"] = False
                localContext["nameSuffix"] = "_WithoutPS"
                TemplatingUtils.renderAsFile("fusedSPTVI_linear_cpp.template", self.m_filename_noPS, localContext)
        else:
            self.m_context["nDof_seq"] = range(0,self.m_context["nDof"])
            self.m_context["gemm_rhs_x"] = gemmNamePad+"_rhs_x"
            self.m_context["gemm_rhs_y"] = gemmNamePad+"_rhs_y"
            self.m_context["gemm_rhs_z"] = gemmNamePad+"_rhs_z"
            self.m_context["gemm_lqi"]   = gemmName+"_lqi"
            self.m_context["gemm_x"] = gemmName+"_lduh_x"
            self.m_context["gemm_y"] = gemmName+"_lduh_y"
            self.m_context["gemm_z"] = gemmName+"_lduh_z"             
            self.m_context["i_seq"] = range(0,self.m_context["nDof"])
            self.m_context["j_seq"] = range(0,self.m_context["nDof"]) if (self.m_context["nDim"] >= 3) else [0]
            
            TemplatingUtils.renderAsFile("fusedSPTVI_nonlinear_cpp.template", self.m_filename, self.m_context)
            
        # generates gemms
        if(self.m_context["useLibxsmm"]):
            self.generateGemms()

    def generateGemms(self):
        l_matmulList = []
        if(self.m_context["isLinear"]):
            if(self.m_context["useFlux"]):
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nVarPad"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"], \
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
                l_matmulList.append(l_matmul)
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nVarPad"] * self.m_context["nDof"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"], \
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
                l_matmulList.append(l_matmul)
                if(self.m_context["nDim"]>=3):
                    l_matmul = MatmulConfig(    # M
                                                self.m_context["nVar"],    \
                                                # N
                                                self.m_context["nDof"],    \
                                                # K
                                                self.m_context["nDof"],    \
                                                # LDA
                                                self.m_context["nVarPad"] * (self.m_context["nDof"]**2), \
                                                # LDB
                                                self.m_context["nDofPad"], \
                                                # LDC
                                                self.m_context["nVarPad"], \
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
                    l_matmulList.append(l_matmul)
            if(self.m_context["useNCP"]):
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nDataPad"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"], \
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
                l_matmulList.append(l_matmul)
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nDataPad"] * self.m_context["nDof"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"] * self.m_context["nDof"], \
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
                l_matmulList.append(l_matmul)
                if(self.m_context["nDim"]>=3):
                    l_matmul = MatmulConfig(    # M
                                                self.m_context["nVar"],    \
                                                # N
                                                self.m_context["nDof"],    \
                                                # K
                                                self.m_context["nDof"],    \
                                                # LDA
                                                self.m_context["nDataPad"] * (self.m_context["nDof"] ** 2), \
                                                # LDB
                                                self.m_context["nDofPad"], \
                                                # LDC
                                                self.m_context["nVarPad"] * (self.m_context["nDof"] ** 2), \
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
                    l_matmulList.append(l_matmul)
        else: #NonLinear
            if(self.m_context["useFlux"]):
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVarPad"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nVarPad"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"], \
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
                l_matmulList.append(l_matmul)
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVarPad"],                             \
                                            # N
                                            self.m_context["nDof"],                             \
                                            # K
                                            self.m_context["nDof"],                             \
                                            # LDA
                                            self.m_context["nVarPad"]* self.m_context["nDof"],     \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"] * self.m_context["nDof"],     \
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
                l_matmulList.append(l_matmul)
                if(self.m_context["nDim"]>=3):
                    l_matmul = MatmulConfig(    # M
                                                self.m_context["nVarPad"],                             \
                                                # N
                                                self.m_context["nDof"],                             \
                                                # K
                                                self.m_context["nDof"],                             \
                                                # LDA
                                                self.m_context["nVarPad"] * (self.m_context["nDof"]**2),     \
                                                # LDB
                                                self.m_context["nDofPad"],                          \
                                                # LDC
                                                self.m_context["nVarPad"] * (self.m_context["nDof"]**2),  \
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
                    l_matmulList.append(l_matmul)
                # (1) MATMUL( lFhi_x(:,:,j,k), TRANSPOSE(Kxi) )
                l_matmul_x = MatmulConfig(  # M
                                            self.m_context["nVar"],       \
                                            # N
                                            self.m_context["nDof"],       \
                                            # K
                                            self.m_context["nDof"],       \
                                            # LDA
                                            self.m_context["nVarPad"],    \
                                            # LDB
                                            self.m_context["nDofPad"],    \
                                            # LDC
                                            self.m_context["nVar"],       \
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
                l_matmulList.append(l_matmul_x)

                # (2) MATMUL( lFhi_y(:,:,i,k), TRANSPOSE(Kxi) )
                l_matmul_y = MatmulConfig(  # M
                                            self.m_context["nVar"],                         \
                                            # N
                                            self.m_context["nDof"],                         \
                                            # K
                                            self.m_context["nDof"],                         \
                                            # LDA
                                            self.m_context["nVarPad"],                      \
                                            # LDB
                                            self.m_context["nDofPad"],                      \
                                            # LDC
                                            self.m_context["nVar"]*self.m_context["nDof"],  \
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
                l_matmulList.append(l_matmul_y)

                if(self.m_context["nDim"]>=3):
                    # (3) MATMUL( lFhi_z(:,:,i,j), TRANSPOSE(Kxi) )
                    l_matmul_z = MatmulConfig(  # M
                                                self.m_context["nVar"],                             \
                                                # N
                                                self.m_context["nDof"],                             \
                                                # K
                                                self.m_context["nDof"],                             \
                                                # LDA
                                                self.m_context["nVarPad"],                          \
                                                # LDB
                                                self.m_context["nDofPad"],                          \
                                                # LDC
                                                self.m_context["nVar"]*(self.m_context["nDof"]**2), \
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
                    l_matmulList.append(l_matmul_z)
            if(self.m_context["useNCP"]):
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nDataPad"] * self.m_context["nDof"], \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"] * self.m_context["nDim"] * self.m_context["nDof"], \
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
                l_matmulList.append(l_matmul)
                l_matmul = MatmulConfig(    # M
                                            self.m_context["nVar"],    \
                                            # N
                                            self.m_context["nDof"],    \
                                            # K
                                            self.m_context["nDof"],    \
                                            # LDA
                                            self.m_context["nDataPad"] * (self.m_context["nDof"] ** 2), \
                                            # LDB
                                            self.m_context["nDofPad"], \
                                            # LDC
                                            self.m_context["nVarPad"] * self.m_context["nDim"] * (self.m_context["nDof"] ** 2), \
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
                l_matmulList.append(l_matmul)
                if(self.m_context["nDim"]>=3):
                    l_matmul = MatmulConfig(    # M
                                                self.m_context["nVar"],    \
                                                # N
                                                self.m_context["nDof"],    \
                                                # K
                                                self.m_context["nDof"],    \
                                                # LDA
                                                self.m_context["nDataPad"] * (self.m_context["nDof"] ** 3), \
                                                # LDB
                                                self.m_context["nDofPad"], \
                                                # LDC
                                                self.m_context["nVarPad"] * self.m_context["nDim"] * (self.m_context["nDof"] ** 3), \
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
                    l_matmulList.append(l_matmul)
            l_matmul = MatmulConfig(    # M
                                        self.m_context["nVar"],                             \
                                        # N
                                        self.m_context["nDof"],                             \
                                        # K
                                        self.m_context["nDof"],                             \
                                        # LDA
                                        self.m_context["nVarPad"]*(self.m_context["nDof"]**self.m_context["nDim"]), \
                                        # LDB
                                        self.m_context["nDofPad"], \
                                        # LDC
                                        self.m_context["nVarPad"], \
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
            l_matmulList.append(l_matmul)
            
        Backend.generateAssemblerCode(self.m_filename_asm, l_matmulList)
