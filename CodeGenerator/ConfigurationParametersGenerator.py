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


import TemplatingUtils


class ConfigurationParametersGenerator:
    m_context = {}

    # name of generated output file
    m_filename = "ConfigurationParameters.cpph"

    
    def __init__(self, i_config):
        self.m_context = i_config


    def generateCode(self):
        self.m_context["isLinearCText"] = "true" if self.m_context["isLinear"] else "false" #c++ true/false instead of True/False

        TemplatingUtils.renderAsFile("configurationParameters_cpph.template", self.m_filename, self.m_context)
