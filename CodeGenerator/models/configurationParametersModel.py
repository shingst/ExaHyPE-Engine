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
# Generate a ccph with getter to the parameters used by the code generator
#


from .abstractModelBaseClass import AbstractModelBaseClass


class ConfigurationParametersModel(AbstractModelBaseClass):

    def generateCode(self):
        self.context["isLinearCText"]          = "true" if self.context["isLinear"]    else "false" #c++ true/false instead of True/False
        self.context["useFluxFactor"]          = 1 if self.context["useFlux"]          else 0
        self.context["useNCPFactor"]           = 1 if self.context["useNCP"]           else 0
        self.context["useSourceFactor"]        = 1 if self.context["useSource"]        else 0
        self.context["usePointSourcesFactor"]  = 1 if self.context["usePointSources"]  else 0
        self.context["useMaterialParamFactor"] = 1 if self.context["useMaterialParam"] else 0
        self.context["useSourceOrNCPFactor"]   = 1 if self.context["useSourceOrNCP"]   else 0
        
        self.render("configurationParameters_cpph.template", "ConfigurationParameters.cpph")
