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
from configuration import Configuration
from controller import Controller

def main():
    """"Check python version, call the controller and run it"""
    # check version. Python 3.3 required
    requiredVersion = (3,3)
    currentVersion  = sys.version_info
    if(requiredVersion > currentVersion):
        sys.exit("Requires Python 3.3 or newer. Abort.")
    
    # generate absolute paths to project and LIBXSMM
    dir = os.path.dirname(__file__)
    absolutePathToRoot =  os.path.abspath(os.path.join(dir,Configuration.pathToExaHyPERoot))
    absolutePathToLibxsmm =  os.path.abspath(os.path.join(dir,Configuration.pathToLibxsmmGemmGenerator))
    
    # create controller and run it (input parsed with argparse)
    control = Controller(absolutePathToRoot, absolutePathToLibxsmm, Configuration.simdWidth)
    control.generateCode()

    
if __name__ == "__main__":
    # execute only if run as a script
    main()
