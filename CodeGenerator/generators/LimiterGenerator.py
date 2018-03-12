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


import Backend
from utils import TemplatingUtils
from utils.MatmulConfig import MatmulConfig


class LimiterGenerator:
    m_context = {}

    # name of generated output file
    m_filename = "limiter.cpp"


    def __init__(self, i_context):
        self.m_context = i_context
        

    def generateCode(self):
        if(not self.m_context['useLimiter']):
            return None
        self.m_context["gemm_dg2fv"] = "gemm_"+str(self.m_context["nVar"])+"_"+str(self.m_context["nDofLim"])+"_"+str(self.m_context["nDof"])+"_dg2fv"
        self.m_context["gemm_fv2dg"] = "gemm_"+str(self.m_context["nVar"]) +"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDofLim"])+"_fv2dg"
        self.m_context["gemm_uh2lob"] = "gemm_"+str(self.m_context["nVar"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])+"_uh2lob"
        
        TemplatingUtils.renderAsFile("limiter_cpp.template", self.m_filename, self.m_context)
        # generates gemms
        if(self.m_context["useLibxsmm"]):
            self.generateGemms()

    def generateGemms(self):
        # define a sequence of matmul configs
        l_matmulList = []

        #-----------------------------
        # implementation file
        #-----------------------------        
        l_dg2fv = MatmulConfig(  # M
                                    self.m_context["nVar"],      \
                                    # N
                                    self.m_context["nDofLim"],       \
                                    # K
                                    self.m_context["nDof"],       \
                                    # LDA
                                    self.m_context["nData"],   \
                                    # LDB
                                    self.m_context["nDofPad"],    \
                                    # LDC
                                    self.m_context["nVar"],   \
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
        l_matmulList.append(l_dg2fv)
        l_fv2dg = MatmulConfig(  # M
                                    self.m_context["nVar"],      \
                                    # N
                                    self.m_context["nDof"],       \
                                    # K
                                    self.m_context["nDofLim"],       \
                                    # LDA
                                    self.m_context["nVar"],   \
                                    # LDB
                                    self.m_context["nDofLimPad"],    \
                                    # LDC
                                    self.m_context["nData"],   \
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
        l_matmulList.append(l_fv2dg)
        l_uh2lob = MatmulConfig(  # M
                                    self.m_context["nVar"],      \
                                    # N
                                    self.m_context["nDof"],       \
                                    # K
                                    self.m_context["nDof"],       \
                                    # LDA
                                    self.m_context["nData"],   \
                                    # LDB
                                    self.m_context["nDofPad"],    \
                                    # LDC
                                    self.m_context["nVar"],   \
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
        l_matmulList.append(l_uh2lob)
        
        
        
        Backend.generateAssemblerCode("asm_"+self.m_filename, l_matmulList)
