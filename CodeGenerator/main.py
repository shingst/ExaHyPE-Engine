#!/bin/env python
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
# Starting point of the code generator
#
# @note
# requires python3

import sys
import os
import controller

def main():
    # --------------------------------------------------------
    # Configuration parameters
    # --------------------------------------------------------
    pathFromHereToExaHyPERoot = "../" #path to the root of ExaHyPe from this file
    pathToLibxsmmGemmGenerator = os.path.join("dependencies", "libxsmm", "bin", "libxsmm_gemm_generator") #path to the gemm generator from this file

    simdWidth = { "noarch" : 1,
                  "wsm"    : 2,
                  "snb"    : 4,
                  "hsw"    : 4,
                  "knc"    : 8,
                  "knl"    : 8 
                }

    # --------------------------------------------------------
    # Require python3.3
    # --------------------------------------------------------
    requiredVersion = (3,3)
    currentVersion  = sys.version_info
    if(requiredVersion > currentVersion):
        sys.exit("Requires Python 3.3 or newer. Abort.")
    
    # Absolute path to project
    dir = os.path.dirname(__file__)
    absolutePathToRoot =  os.path.abspath(os.path.join(dir,pathFromHereToExaHyPERoot))
    absolutePathToLibxsmm =  os.path.abspath(os.path.join(dir,pathToLibxsmmGemmGenerator))
    
    # create controller and run it (input parsed with argparse)
    control = controller.Controller(absolutePathToRoot, absolutePathToLibxsmm, simdWidth)
    control.generateCode()

    
if __name__ == "__main__":
    # execute only if run as a script
    main()
