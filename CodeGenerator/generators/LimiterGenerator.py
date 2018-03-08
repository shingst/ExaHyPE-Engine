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


from utils import TemplatingUtils


class LimiterGenerator:
    m_context = {}

    # name of generated output file
    m_filename = "limiter.cpp"


    def __init__(self, i_context):
        self.m_context = i_context
        

    def generateCode(self):  
        TemplatingUtils.renderAsFile("limiter_cpp.template", self.m_filename, self.m_context)
