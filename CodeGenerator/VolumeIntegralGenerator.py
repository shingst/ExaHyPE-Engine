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
# Generates the VolumeIntegral kernel
#


import Backend
import TemplatingUtils
from MatmulConfig import MatmulConfig


class VolumeIntegralGenerator:
    m_context = {}
    
    # name of generated output file
    m_filename = "volumeIntegral.cpp"


    def __init__(self, i_context):
        self.m_context = i_context


    def generateCode(self):
        if(self.m_context["isLinear"]):
            TemplatingUtils.renderAsFile("volumeIntegralLinear_cpp.template", self.m_filename, self.m_context)
        else:
            # initialize context
            gemmName = "gemm_"+str(self.m_context["nVar"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])
            self.m_context["gemm_x"] = gemmName+"_lduh_x"
            self.m_context["gemm_y"] = gemmName+"_lduh_y"
            self.m_context["gemm_z"] = gemmName+"_lduh_z"             
            self.m_context["i_seq"] = range(0,self.m_context["nDof"])
            self.m_context["j_seq"] = range(0,self.m_context["nDof"]) if (self.m_context["nDim"] >= 3) else [0]
            
            # render template
            if(self.m_context["noTimeAveraging"]):
                TemplatingUtils.renderAsFile("volumeIntegralNonLinear_noTimeAveraging_cpp.template", self.m_filename, self.m_context)
            else:
                TemplatingUtils.renderAsFile("volumeIntegralNonLinear_cpp.template", self.m_filename, self.m_context)
            
            # generates gemms
            if(self.m_context["useLibxsmm"] and self.m_context["useFlux"]): #no gemm if no flux
                self.generateNonlinearGemms() 


    def generateNonlinearGemms(self):
        # Compute
        # (1) lduh(:,:,j,k) = lduh(:,:,j,k) + MATMUL( lFhi_x(:,:,j,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(1)
        # (2) lduh(:,i,:,k) = lduh(:,i,:,k) + MATMUL( lFhi_y(:,:,i,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(2)
        # (3) lduh(:,i,j,:) = lduh(:,i,j,:) + MATMUL( lFhi_z(:,:,i,j), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(3)

        # define a sequence of matmul configs
        l_matmulList = []

        #-----------------------------
        # implementation file
        #-----------------------------

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

        Backend.generateAssemblerCode("asm_"+self.m_filename, l_matmulList)
