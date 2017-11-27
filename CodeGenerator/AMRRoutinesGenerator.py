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


import Backend
import TemplatingUtils
from MatmulConfig import MatmulConfig


class AMRRoutinesGenerator:
    m_context = {}
    
    # name of generated output file
    m_filename = "amrRoutines.cpp"


    def __init__(self, i_context):
        self.m_context = i_context


    def generateCode(self):
        self.m_context["gemm_face_Q"] = "gemm_"+str(self.m_context["nData"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])+"_face_Q"
        self.m_context["gemm_face_F"] = "gemm_"+str(self.m_context["nVar"]) +"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])+"_face_F"
        self.m_context["gemm_volume"] = "gemm_"+str(self.m_context["nData"])+"_"+str(self.m_context["nDof"])+"_"+str(self.m_context["nDof"])+"_volume"
        
        TemplatingUtils.renderAsFile("amrRoutines_cpp.template", self.m_filename, self.m_context)
        # generates gemms
        if(self.m_context["useLibxsmm"]):
            self.generateGemms()
    
    def generateGemms(self):
        # define a sequence of matmul configs
        l_matmulList = []

        #-----------------------------
        # implementation file
        #-----------------------------        
        l_face_Q = MatmulConfig(  # M
                                    self.m_context["nData"],      \
                                    # N
                                    self.m_context["nDof"],       \
                                    # K
                                    self.m_context["nDof"],       \
                                    # LDA
                                    self.m_context["nDataPad"],   \
                                    # LDB
                                    self.m_context["nDofPad"],    \
                                    # LDC
                                    self.m_context["nDataPad"],   \
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
        l_matmulList.append(l_face_Q)
        l_face_F = MatmulConfig(  # M
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
                                    self.m_context["nVarPad"],    \
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
        l_matmulList.append(l_face_F)
        l_volume = MatmulConfig(  # M
                                    self.m_context["nData"],       \
                                    # N
                                    self.m_context["nDof"],       \
                                    # K
                                    self.m_context["nDof"],       \
                                    # LDA
                                    self.m_context["nData"],      \
                                    # LDB
                                    self.m_context["nDofPad"],    \
                                    # LDC
                                    self.m_context["nData"],      \
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
        l_matmulList.append(l_volume)
        
        Backend.generateAssemblerCode("asm_"+self.m_filename, l_matmulList)
