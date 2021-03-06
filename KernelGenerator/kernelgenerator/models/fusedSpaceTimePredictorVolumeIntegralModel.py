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

from ..utils import MatmulConfig


class FusedSpaceTimePredictorVolumeIntegralModel(AbstractModelBaseClass):
    
    def generateCode(self):
        self.context["nVarMinusOne_seq"] = range(self.context["nVar"] - 1)
        self.context["nDataMinusOne_seq"] = range(self.context["nData"] - 1)
        
        if self.context["isLinear"]:

            if self.context["useSplitCK"]:
                if self.context["useVectPDE"]:
                    if self.context["useAoSoA2"]:
                        if self.context["singlePrecisionSTP"]:
                            template = "SP_fusedSPTVI_linear_split_ck_aosoa2_cpp.template"
                        else:
                            template = "fusedSPTVI_linear_split_ck_aosoa2_cpp.template"
                    else:
                        template = "fusedSPTVI_linear_split_ck_vect_cpp.template"
                else:
                    template = "fusedSPTVI_linear_split_ck_cpp.template"

                self.render(("aderdg", template), "fusedSpaceTimePredictorVolumeIntegral.cpp")
                
                if self.context["usePointSources"]:
                    self.context["usePointSources"] = False
                    self.context["nameSuffix"] = "_WithoutPS"
                    
                    self.render(("aderdg", template), "fusedSpaceTimePredictorVolumeIntegral_WithoutPS.cpp")
            
            else:
                # size of the tmpArray
                self.context["tmpArraySize"] = max((self.context["nDof"]*self.context["nVarPad"] if self.context["useFlux"]          else 0), \
                                                   (self.context["nDim"]*self.context["nVarPad"] if self.context["useNCP"]           else 0))
                
                self.render(("aderdg", "fusedSPTVI_linear_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral.cpp")
                
                if self.context["usePointSources"]:
                    self.context["usePointSources"] = False
                    self.context["nameSuffix"] = "_WithoutPS"
                    
                    self.render(("aderdg", "fusedSPTVI_linear_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral_WithoutPS.cpp")
                
        else: #nonlinear
            self.context["nDof_seq"] = range(0,self.context["nDof"])
            self.context["i_seq"] = range(0,self.context["nDof"])
            self.context["j_seq"] = range(0,self.context["nDof"]) if (self.context["nDim"] >= 3) else [0]
            if self.context["predictorRecompute"]:
                if self.context["useVectPDE"]:
                    if self.context["useAoSoA2"]:
                        self.render(("aderdg", "fusedSPTVI_nonlinear_mem_aosoa2_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral.cpp")
                    else:
                        self.render(("aderdg", "fusedSPTVI_nonlinear_mem_vect_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral.cpp")
                else:
                    self.render(("aderdg", "fusedSPTVI_nonlinear_mem_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral.cpp")
            else:
                self.render(("aderdg", "fusedSPTVI_nonlinear_cpp.template"), "fusedSpaceTimePredictorVolumeIntegral.cpp")
       
        # return required gemms
        return {"matmulList": self.context["matmulConfigs"].values(), "fileName": "asm_fstpvi.c"}
    
    
    def buildGemmsConfig(self):
        # MatmulConfig: M, N, K, LDA, LDB, LDC, alpha, beta, Align A, Align B, Align C, name, prefetching, type
        self.context["matmulConfigs"] = {}
        # shortcut
        nVar     = self.context["nVar"]
        nVarPad  = self.context["nVarPad"]
        nDataPad = self.context["nDataPad"]
        nData    = self.context["nData"]
        nDof     = self.context["nDof"]
        nDof2    = nDof*nDof
        nDof3    = nDof2*nDof
        nDof3D   = self.context["nDof3D"]
        nDofPad  = self.context["nDofPad"]
        nDof2Pad = self.context["nDof2Pad"]
        nDim     = self.context["nDim"]
        
        prec = "DP" if not self.context["useSinglePrecision"] else "SP" # only used for Linear SplitCK AoSoA2
        
    # LINEAR
        
        if self.context["isLinear"]:
            if self.context["useSplitCK"]:
                if self.context["useVectPDE"]: 
                    if self.context["useAoSoA2"]: #split_ck aosoa2
                    
        # Linear SplitCK AoSoA2
                    
                        if self.context["useFlux"]:
                            if self.context["useMaterialParam"]:
                                self.context["matmulConfigs"]["flux_x_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDofPad , nDof, nDof    , 1, 0, 1, 1, 1, "flux_x_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                                self.context["matmulConfigs"]["flux_y_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDof    , nDofPad, nDof , 1, 0, 1, 1, 1, "flux_y_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                                if self.context["nDim"]>=3:
                                    self.context["matmulConfigs"]["flux_z_sck_aosoa2"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar , nDofPad, nDof2Pad*nVar , 1, 0, 1, 1, 1, "flux_z_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                            else:
                                self.context["matmulConfigs"]["flux_x_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDofPad , nDof, nDof    , 1, 1, 1, 1, 1, "flux_x_sck_aosoa2", precision=prec)
                                self.context["matmulConfigs"]["flux_y_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDof    , nDofPad, nDof , 1, 1, 1, 1, 1, "flux_y_sck_aosoa2", precision=prec)
                                if self.context["nDim"]>=3:
                                    self.context["matmulConfigs"]["flux_z_sck_aosoa2"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar , nDofPad, nDof2Pad*nVar , 1, 1, 1, 1, 1, "flux_z_sck_aosoa2", precision=prec)
                        self.context["matmulConfigs"]["gradQ_x_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDofPad , nDof, nDof               , 1, 0, 1, 1, 1, "gradQ_x_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                        self.context["matmulConfigs"]["gradQ_y_sck_aosoa2"] =     MatmulConfig(nDof, nDof, nDof, nDof , nDofPad, nDof               , 1, 0, 1, 1, 1, "gradQ_y_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["gradQ_z_sck_aosoa2"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar, nDofPad, nDof2Pad*nVar, 1, 0, 1, 1, 1, "gradQ_z_sck_aosoa2", precision=prec) # beta, 0 => overwrite C
                    
        # Linear SplitCK vect (AoSoA)
                    
                    else:
                        if self.context["useFlux"]:
                            if self.context["useMaterialParam"]:
                                self.context["matmulConfigs"]["flux_x_sck_vect"] =       MatmulConfig(nDofPad, nVar, nDof, nDofPad , nDofPad, nDofPad      , 1, 0, 1, 1, 1, "flux_x_sck_vect") # beta, 0 => overwrite C
                                self.context["matmulConfigs"]["flux_y_or_z_sck_vect"] =       MatmulConfig(nDofPad*nVar, nDof, nDof, nDofPad*nVar , nDofPad, nDofPad*nVar      , 1, 0, 1, 1, 1, "flux_y_or_z_sck_vect") # no moving slices, beta, 0 => overwrite C
                            else:
                                self.context["matmulConfigs"]["flux_x_sck_vect"] =  MatmulConfig(nDofPad, nVar, nDof, nDofPad      , nDofPad, nDofPad      , 1, 1, 1, 1, 1, "flux_x_sck_vect", prefetchOutput="C")
                                self.context["matmulConfigs"]["flux_y_sck_vect"] =       MatmulConfig(nDofPad*nVar, nDof, nDof, nDofPad*nVar , nDofPad, nDofPad*nVar , 1, 1, 1, 1, 1, "flux_y_sck_vect", prefetchOutput="C")
                                if self.context["nDim"]>=3:
                                    self.context["matmulConfigs"]["flux_z_sck_vect"] =   MatmulConfig(nDofPad*nVar, nDof, nDof, nDofPad*nVar , nDofPad, nDofPad*nVar*nDof, 1, 1, 1, 1, 1, "flux_z_sck_vect", prefetchOutput="C")
                        self.context["matmulConfigs"]["gradQ_x_sck_vect"] =     MatmulConfig(nDofPad, nVar*nDof*nDof3D, nDof, nDofPad , nDofPad, nDofPad    , 1, 0, 1, 1, 1, "gradQ_x_sck_vect") # no moving slices, beta, 0 => overwrite C
                        self.context["matmulConfigs"]["gradQ_y_sck_vect"] =     MatmulConfig(nDofPad*nVar, nDof, nDof, nDofPad*nVar , nDofPad, nDofPad*nVar , 1, 0, 1, 1, 1, "gradQ_y_sck_vect", prefetchInput="A", prefetchOutput="C") # beta, 0 => overwrite C
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["gradQ_z_sck_vect"] = MatmulConfig(nDofPad*nVar*nDof, nDof, nDof, nDofPad*nVar*nDof, nDofPad, nDofPad*nVar*nDof, 1, 0, 1, 1, 1, "gradQ_z_sck_vect") # no moving slices, beta, 0 => overwrite C
                            
        # Linear SplitCK scalar
                            
                else: 
                    if self.context["useFlux"]:
                        self.context["matmulConfigs"]["flux_x_sck"] =       MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "flux_x_sck")
                        self.context["matmulConfigs"]["flux_y_sck"] =       MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "flux_y_sck")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["flux_z_sck"] =   MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "flux_z_sck")
                    self.context["matmulConfigs"]["gradQ_x_sck"] =          MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad      , 1, 0, 1, 1, 1, "gradQ_x_sck") # beta, 0 => overwrite C
                    self.context["matmulConfigs"]["gradQ_y_sck"] =          MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof , nDofPad, nVarPad*nDof , 1, 0, 1, 1, 1, "gradQ_y_sck") # beta, 0 => overwrite C
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["gradQ_z_sck"] =      MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 0, 1, 1, 1, "gradQ_z_sck", "nopf", "gemm") # beta, 0 => overwrite C
            
        # Linear default
            
            else: 
                if self.context["useFlux"]:
                    self.context["matmulConfigs"]["flux_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad, 1, 0, 1, 1, 1, "flux_x") # beta, 0 => overwrite C
                    self.context["matmulConfigs"]["flux_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof , nDofPad, nVarPad, 1, 0, 1, 1, 1, "flux_y") # beta, 0 => overwrite C
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["flux_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad, 1, 0, 1, 1, 1, "flux_z") # beta, 0 => overwrite C
                if self.context["useNCP"]:
                    self.context["matmulConfigs"]["gradQ_x"] =     MatmulConfig(nVar, nDof, nDof, nDataPad      , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "gradQ_x")
                    self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nVar, nDof, nDof, nDataPad*nDof , nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "gradQ_y")
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["gradQ_z"] = MatmulConfig(nVar, nDof, nDof, nDataPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "gradQ_z")
        
    # NONLINEAR
        
        else:
            if self.context["predictorRecompute"]: # TODO JMG matmuls for gradQ, rhs and lduh are exactly the same...
                if self.context["useVectPDE"]:
                    if self.context["useAoSoA2"]: # aosoa2
                    
        # Nonlinear PredictorRecompute AoSoA2
        
                        if self.context["useFlux"]:
                            self.context["matmulConfigs"]["rhs_x"] =     MatmulConfig(nDof, nDof, nDof, nDofPad, nDof, nDof, 1, 1, 1, 1, 1, "rhs_x") #, prefetchInput="B", prefetchOutput="C")
                            self.context["matmulConfigs"]["rhs_y"] =     MatmulConfig(nDof, nDof, nDof, nDof, nDofPad, nDof, 1, 1, 1, 1, 1, "rhs_y") #, prefetchInput="A", prefetchOutput="C")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["rhs_z"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar, nDofPad, nDof2Pad*nVar, 1, 1, 1, 1, 1, "rhs_z")
                            self.context["matmulConfigs"]["lduh_x"] =     MatmulConfig(nDof, nDof, nDof, nDofPad, nDof, nDof, 1, 1, 1, 1, 1, "lduh_x") #, prefetchInput="B", prefetchOutput="C")
                            self.context["matmulConfigs"]["lduh_y"] =     MatmulConfig(nDof, nDof, nDof, nDof, nDofPad, nDof, 1, 1, 1, 1, 1, "lduh_y") #, prefetchInput="A", prefetchOutput="C")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["lduh_z"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar, nDofPad, nDof2Pad*nVar, 1, 1, 1, 1, 1, "lduh_z")
                        if self.context["useNCP"] or self.context['useViscousFlux']:
                            self.context["matmulConfigs"]["gradQ_x"] =     MatmulConfig(nDof, nDof, nDof, nDofPad, nDof, nDof, 1, 1, 1, 1, 1, "gradQ_x")
                            self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nDof, nDof, nDof, nDof, nDofPad, nDof, 1, 1, 1, 1, 1, "gradQ_y")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["gradQ_z"] = MatmulConfig(nDof2Pad*nVar, nDof, nDof, nDof2Pad*nVar, nDofPad, nDof2Pad*nVar, 1, 1, 1, 1, 1, "gradQ_z")
                        self.context["matmulConfigs"]["lqi"] = MatmulConfig(nDof2Pad, nDof, nDof, nDof2Pad*nVar*nDof3D, nDofPad, nDof2Pad, 1, 0, 1, 1, 1, "lqi") # beta, 0 => overwrite C
                    
        # Nonlinear PredictorRecompute vect (AoSoA)
                    
                    else:
                        if self.context["useFlux"]:
                            self.context["matmulConfigs"]["rhs_x"] =     MatmulConfig(nDofPad,           nVar, nDof, nDofPad,           nDofPad, nDofPad          , 1, 1, 1, 1, 1, "rhs_x", prefetchInput="B", prefetchOutput="C")
                            self.context["matmulConfigs"]["rhs_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "rhs_y", prefetchInput="A", prefetchOutput="C")
                            #self.context["matmulConfigs"]["rhs_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "rhs_y")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["rhs_z"] = MatmulConfig(nDofPad*nVar*nDof, nDof, nDof, nDofPad*nVar*nDof, nDofPad, nDofPad*nVar*nDof , 1, 1, 1, 1, 1, "rhs_z")
                            self.context["matmulConfigs"]["lduh_x"] =     MatmulConfig(nDofPad,           nVar, nDof, nDofPad,           nDofPad, nDofPad          , 1, 1, 1, 1, 1, "lduh_x", prefetchInput="B", prefetchOutput="C")
                            self.context["matmulConfigs"]["lduh_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "lduh_y", prefetchInput="A", prefetchOutput="C")
                            #self.context["matmulConfigs"]["lduh_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "lduh_y")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["lduh_z"] = MatmulConfig(nDofPad*nVar*nDof, nDof, nDof, nDofPad*nVar*nDof, nDofPad, nDofPad*nVar*nDof, 1, 1, 1, 1, 1, "lduh_z")
                        if self.context["useNCP"] or self.context['useViscousFlux']:
                            self.context["matmulConfigs"]["gradQ_x"] =     MatmulConfig(nDofPad,           nVar, nDof, nDofPad,           nDofPad, nDofPad          , 1, 1, 1, 1, 1, "gradQ_x", prefetchInput="B", prefetchOutput="C")
                            self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "gradQ_y", prefetchInput="A", prefetchOutput="C")
                            #self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nDofPad*nVar,      nDof, nDof, nDofPad*nVar,      nDofPad, nDofPad*nVar     , 1, 1, 1, 1, 1, "gradQ_y")
                            if self.context["nDim"]>=3:
                                self.context["matmulConfigs"]["gradQ_z"] = MatmulConfig(nDofPad*nVar*nDof, nDof, nDof, nDofPad*nVar*nDof, nDofPad, nDofPad*nVar*nDof, 1, 1, 1, 1, 1, "gradQ_z")
                        self.context["matmulConfigs"]["lqi"] = MatmulConfig(nDofPad*nVar, nDof, nDof, nDofPad*nVar*nDof*nDof3D, nDofPad, nDofPad*nVar, 1, 0, 1, 1, 1, "lqi") # beta, 0 => overwrite C
                
        # Nonlinear PredictorRecompute scalar
                
                else:
                    if self.context["useFlux"]:
                        self.context["matmulConfigs"]["rhs_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "rhs_x")
                        self.context["matmulConfigs"]["rhs_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof , nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "rhs_y")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["rhs_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "rhs_z")
                        self.context["matmulConfigs"]["lduh_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad, nDofPad, nVarPad            , 1, 1, 1, 1, 1, "lduh_x")
                        self.context["matmulConfigs"]["lduh_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof, nDofPad, nVarPad*nDof  , 1, 1, 1, 1, 1, "lduh_y")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["lduh_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "lduh_z")
                    if self.context["useNCP"] or self.context['useViscousFlux']:
                        self.context["matmulConfigs"]["gradQ_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad , nDofPad, nVarPad           , 1, 1, 1, 1, 1, "gradQ_x")
                        self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof, nDofPad, nVarPad*nDof  , 1, 1, 1, 1, 1, "gradQ_y")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["gradQ_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "gradQ_z")
                    self.context["matmulConfigs"]["lqi"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*(nDof**nDim), nDofPad, nVarPad, 1, 0, 1, 1, 1, "lqi") # beta, 0 => overwrite C
            
        # Nonlinear default
            
            else:
                if self.context["useFlux"]:
                    self.context["matmulConfigs"]["rhs_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad      , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "rhs_x")
                    self.context["matmulConfigs"]["rhs_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof , nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "rhs_y")
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["rhs_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "rhs_z")

                    self.context["matmulConfigs"]["lduh_x"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad, nDofPad, nVarPad      , 1, 1, 1, 1, 1, "lduh_x")
                    self.context["matmulConfigs"]["lduh_y"] =     MatmulConfig(nVarPad, nDof, nDof, nVarPad, nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "lduh_y")
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["lduh_z"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "lduh_z")
                    if self.context["useCERKGuess"]:
                        self.context["matmulConfigs"]["gradF_x_RKLoop"] =     MatmulConfig(nVar, nDof, nDof, nVarPad      , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "gradF_x_RKLoop")
                        self.context["matmulConfigs"]["gradF_y_RKLoop"] =     MatmulConfig(nVar, nDof, nDof, nVarPad*nDof , nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "gradF_y_RKLoop")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["gradF_z_RKLoop"] = MatmulConfig(nVar, nDof, nDof, nVarPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "gradF_z_RKLoop")
                if self.context["useNCP"] or self.context['useViscousFlux']:
                    self.context["matmulConfigs"]["gradQ_x"] =     MatmulConfig(nVar, nDof, nDof, nDataPad , nDofPad, nVarPad      , 1, 1, 1, 1, 1, "gradQ_x")
                    self.context["matmulConfigs"]["gradQ_y"] =     MatmulConfig(nVar, nDof, nDof, nDataPad*nDof, nDofPad, nVarPad*nDof , 1, 1, 1, 1, 1, "gradQ_y")
                    if self.context["nDim"]>=3:
                        self.context["matmulConfigs"]["gradQ_z"] = MatmulConfig(nVar, nDof, nDof, nDataPad*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 1, 1, 1, "gradQ_z")
                    if self.context["useCERKGuess"]:
                        self.context["matmulConfigs"]["gradQ_x_RKLoop"] =     MatmulConfig(nVar, nDof, nDof, nData      , nDofPad, nVarPad      , 1, 1, 0, 1, 1, "gradQ_x_RKLoop")
                        self.context["matmulConfigs"]["gradQ_y_RKLoop"] =     MatmulConfig(nVar, nDof, nDof, nData*nDof , nDofPad, nVarPad*nDof , 1, 1, 0, 1, 1, "gradQ_y_RKLoop")
                        if self.context["nDim"]>=3:
                            self.context["matmulConfigs"]["gradQ_z_RKLoop"] = MatmulConfig(nVar, nDof, nDof, nData*nDof2, nDofPad, nVarPad*nDof2, 1, 1, 0, 1, 1, "gradQ_z_RKLoop")
                self.context["matmulConfigs"]["lqi"] = MatmulConfig(nVarPad, nDof, nDof, nVarPad*(nDof**nDim), nDofPad, nVarPad, 1, 0, 1, 1, 1, "lqi") # beta, 0 => overwrite C
