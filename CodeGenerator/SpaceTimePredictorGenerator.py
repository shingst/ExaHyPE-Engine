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


import Backend
import TemplatingUtils
from MatmulConfig import MatmulConfig


class SpaceTimePredictorGenerator:
    m_context = {}

    # name of generated output file
    m_filename_picard       = "picard.cpp"
    m_filename_predictor    = "predictor.cpp"
    m_filename_extrapolator = "extrapolatedPredictor.cpp"
    m_filename_linear       = "spaceTimePredictorLinear.cpp"
    
    m_filename_asm_picard   = "asm_picard" 

    
    def __init__(self, i_config):
        self.m_context = i_config


    def generateCode(self):
        if(self.m_context["isLinear"]):
            TemplatingUtils.renderAsFile("spaceTimePredictorLinear_cpp.template", self.m_filename_linear, self.m_context)
        else:
            self.m_context["nDof_seq"] = range(0,self.m_context["nDof"])
            gemmName = "gemm_"+str(self.m_context["nVar"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])
            self.m_context["gemm_rhs_x"] = gemmName+"_rhs_x"
            self.m_context["gemm_rhs_y"] = gemmName+"_rhs_y"
            self.m_context["gemm_rhs_z"] = gemmName+"_rhs_z"
            self.m_context["gemm_gradQ_x"] = gemmName+"_gradQ_x"
            self.m_context["gemm_gradQ_y"] = gemmName+"_gradQ_y"
            self.m_context["gemm_gradQ_z"] = gemmName+"_gradQ_z"
            self.m_context["gemm_lqi"]   = gemmName+"_lqi"
            
            TemplatingUtils.renderAsFile("spaceTimePredictor_picard_cpp.template", self.m_filename_picard, self.m_context)
            if(self.m_context["noTimeAveraging"]):  
                TemplatingUtils.renderAsFile("spaceTimePredictor_extrapolator_noTimeAveraging_cpp.template", self.m_filename_extrapolator, self.m_context)
            else:
                TemplatingUtils.renderAsFile("spaceTimePredictor_predictor_cpp.template", self.m_filename_predictor, self.m_context)
                TemplatingUtils.renderAsFile("spaceTimePredictor_extrapolator_cpp.template", self.m_filename_extrapolator, self.m_context)
            
            # generates gemms
            if(self.m_context["useLibxsmm"]):
                self.generateGemms()

    def generateGemms(self):
        l_matmulList = []
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
                                        self.m_context["nVar"],                             \
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
                                            self.m_context["nVar"],                             \
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
        Backend.generateAssemblerCode(self.m_filename_asm_picard, l_matmulList)
