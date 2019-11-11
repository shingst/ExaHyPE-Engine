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
# Generates the RiemannSolver Kernel
#
# Call user function flux, ncp, eigenvalue


from .abstractModelBaseClass import AbstractModelBaseClass
import math


class RiemannModel(AbstractModelBaseClass):

    def generateCode(self):
        self.context["sqrt_half_Pi"] = math.sqrt(math.pi/2)
        
        if self.context["kernelType"] == "aderdg":
            if(self.context["isLinear"]):
                self.render(("aderdg", "riemannSolverLinear_cpp.template"),    "riemannSolver.cpp")
            else:
                self.render(("aderdg", "riemannSolverNonLinear_cpp.template"), "riemannSolver.cpp")
        elif self.context["kernelType"] == "fv":
            self.render(("fv", "riemannSolver_cpp.template"), "riemannSolver.cpp")
