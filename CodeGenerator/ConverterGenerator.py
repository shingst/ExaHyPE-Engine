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


import TemplatingUtils


class ConverterGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "converter"


    def __init__(self, i_context):
        self.m_context = i_context
        

    def generateCode(self):
        self.m_context["noVarPadding"] = self.m_context["nVarPad"] == self.m_context["nVar"]
        
        TemplatingUtils.renderAsFile("converter_h.template",   self.m_filenameRoot+".h",   self.m_context)
        TemplatingUtils.renderAsFile("converter_cpp.template", self.m_filenameRoot+".cpp", self.m_context)
