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
# Generate the Quadrature matrices (nodes + weights) used by the solver
# Use pure python matrices operations from Utils
#


import Backend
from utils import TemplatingUtils
from utils import MathsUtils #matrix operation and build functions


class QuadratureGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "Quadrature"
    
    # quadrature nodes and weights mapped onto [0,1]
    m_QuadratureWeights       = []
    

    def __init__(self, i_context):
        self.m_context = i_context
        self.m_QuadratureWeights, self.m_QuadratureNodes = MathsUtils.getGaussLegendre(self.m_context["nDof"])
        self.m_OtherQuadratureWeights, self.m_OtherQuadratureNodes = MathsUtils.getGaussLobatto(self.m_context["nDof"])
                

    def generateCode(self):
        self.m_context["quadratureType"] = "Gauss-Legendre"
        l_weightsVector = MathsUtils.vectorPad(self.m_QuadratureWeights, self.m_context["nDofPad"] - self.m_context["nDof"])
        self.m_context["weights1"] = l_weightsVector
        self.m_context["w1Size"] = len(self.m_context["weights1"])
        self.m_context["w1_seq"] = range(self.m_context["w1Size"])
        if(self.m_context["nDim"] == 2):
            # weightsVector is QuadratureWeights itself
            l_weightsVector      = MathsUtils.vectorPad(self.m_QuadratureWeights, Backend.getPadWidth(len(self.m_QuadratureWeights)))
            self.m_context["weights2"] = l_weightsVector
            self.m_context["w2Size"] = len(self.m_context["weights2"])
            self.m_context["w2_seq"] = range(self.m_context["w2Size"])

            # all combinations of two weights, written as an 1D array
            l_weightsVector = [self.m_QuadratureWeights[i] * self.m_QuadratureWeights[j] for i in range(self.m_context["nDof"]) for j in range(self.m_context["nDof"])]
            l_weightsVector = MathsUtils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context["weights3"] = l_weightsVector
            self.m_context["w3Size"] = len(self.m_context["weights3"])
            self.m_context["w3_seq"] = range(self.m_context["w3Size"])

        elif(self.m_context["nDim"] == 3):
            # all combinations of two weights, written as an 1D array
            l_weightsVector = [self.m_QuadratureWeights[i] * self.m_QuadratureWeights[j] for i in range(self.m_context["nDof"]) for j in range(self.m_context["nDof"])]
            l_weightsVector = MathsUtils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context["weights2"] = l_weightsVector
            self.m_context["w2Size"] = len(self.m_context["weights2"])
            self.m_context["w2_seq"] = range(self.m_context["w2Size"])

            # all combination of three weights, written as an 1D array
            l_weightsVector = [self.m_QuadratureWeights[i] * self.m_QuadratureWeights[j] * self.m_QuadratureWeights[k] for i in range(self.m_context["nDof"]) for j in range(self.m_context["nDof"]) for k in range(self.m_context["nDof"])]
            l_weightsVector = MathsUtils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context["weights3"] = l_weightsVector
            self.m_context["w3Size"] = len(self.m_context["weights3"])
            self.m_context["w3_seq"] = range(self.m_context["w3Size"])
            
        else:
            print("WeightsGenerator.__generateWeightsCombinations(): nDim not supported")

        self.m_context["QuadratureWeights"], self.m_context["xGPN"] = MathsUtils.getGaussLegendre(self.m_context["nDof"])
        self.m_context["GPN_seq"] = range(self.m_context["nDof"])
        
        if(self.m_context["useLimiter"]):
            l_padSize = self.m_context["nDofPad"] - self.m_context["nDof"]
            uh2lob = MathsUtils.assembleQuadratureConversion(self.m_QuadratureNodes, self.m_OtherQuadratureNodes, self.m_context["nDof"]) #TODO JMG adapt when allowing Lobatto as node
            self.m_context["uh2lob"] = MathsUtils.matrixPadAndFlatten_ColMajor(uh2lob,l_padSize)
            self.m_context["uh2lobSize"] = len(self.m_context["uh2lob"])
            self.m_context["uh2lob_seq"] = range(self.m_context["uh2lobSize"])
            
            l_padSize = self.m_context["nDofPad"] - self.m_context["nDof"]
            dg2fv = MathsUtils.assembleDGToFV(self.m_QuadratureNodes, self.m_QuadratureWeights, self.m_context["nDof"], self.m_context["nDofLim"])
            self.m_context["dg2fv"] = MathsUtils.matrixPadAndFlatten_ColMajor(dg2fv,l_padSize)
            self.m_context["dg2fvSize"] = len(self.m_context["dg2fv"])
            self.m_context["dg2fv_seq"] = range(self.m_context["dg2fvSize"])
            
            l_padSize = self.m_context["nDofLimPad"] - self.m_context["nDofLim"]
            fv2dg = MathsUtils.assembleFVToDG(dg2fv, self.m_QuadratureWeights, self.m_context["nDof"], self.m_context["nDofLim"])
            self.m_context["fv2dg"] = MathsUtils.matrixPadAndFlatten_ColMajor(fv2dg,l_padSize)
            self.m_context["fv2dgSize"] = len(self.m_context["fv2dg"])
            self.m_context["fv2dg_seq"] = range(self.m_context["fv2dgSize"])
            
        
        #generate files 
        TemplatingUtils.renderAsFile("Quadrature_h.template",   self.m_filenameRoot+".h",   self.m_context)
        TemplatingUtils.renderAsFile("Quadrature_cpp.template", self.m_filenameRoot+".cpp", self.m_context)
