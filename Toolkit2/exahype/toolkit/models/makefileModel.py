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
# Generates the kernel for the mapping of the DG polynomial
# onto the [0,1] and calls the user-defined function.
#
# Call the user function adjustPointSolution
#

from .abstractModelBaseClass import AbstractModelBaseClass

import os

# import Configuration for the alignements
import sys
sys.path.append("..")
from configuration import Configuration


class MakefileModel(AbstractModelBaseClass):

    def setupContext(self):
        """Setup the local context.
        Required to ensure both generateCode and getOutputMessage have a properly set context
        """
        if "setupContextFlag" in self.context:
            return
        
        # dependencies
        self.context["useSharedMem"]      = "shared_memory" in self.context;
        self.context["useDistributedMem"] = "distributed_memory" in self.context;
        if not "TBB_INC" in os.environ:
            print("WARNING: environment variable TBB_INC not set but required if code is built with TBB");
        if not "TBB_SHLIB" in os.environ:
            print("WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB");
        self.context["useIpcm"]   = False # TODO
        self.context["useLikwid"] = False # TODO
        self.context["likwidInc"] = ""    # TODO
        
        # architecture alignement
        self.context["alignment"] = Configuration.alignmentPerArchitectures[self.context["architecture"]]
        
        # kernels
        useOptKernel = False
        useFortran   = False
        for solver in self.context["solvers"]:
            if "aderdg_kernel" in solver:
                useOptKernel = useOptKernel or solver["aderdg_kernel"].get("type","generic")=="optimised"
                useFortran   = useFortran or solver["aderdg_kernel"].get("language","C")=="Fortran"
            if "fv_kernel" in solver:
                useOptKernel = useOptKernel or solver["fv_kernel"].get("type","generic")=="optimised"
                useFortran   = useFortran or solver["fv_kernel"].get("language","C")=="Fortran"
        self.context["useOptKernel"] = useOptKernel
        self.context["useFortran"]   = useFortran
        
        self.context["setupContextFlag"] = True #set a flag to not run twice

    def generateCode(self):
        self.setupContext()
        return self.render("Makefile.template", "Makefile") #return path to generated file
    
    
    def getOutputMessage(self):
        self.setupContext()
        return self.renderAsString("buildMessage.template")