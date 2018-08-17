import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required

    #TODO Add missing path
    # absolute path to ExaHyPE's root (we need absolute paths in generated Makefile)
    # TODO: This is flawed, we have some paths in the specfile but not this one. Why?
    pathToExaHyPERoot          = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    
    # absolute path to jinja2
    pathToJinja2               = os.path.abspath(os.path.join(os.path.dirname(__file__),"..", "..", "dependencies", "jinja"))
    
    # absolute path to markupsafe
    pathToMarkupsafe           = os.path.abspath(os.path.join(os.path.dirname(__file__),"..", "..", "dependencies", "markupsafe"))
    
    alignmentPerArchitectures  = {
        "noarch" : 16,
        "snb"    : 32, 
        "hsw"    : 32, 
        "knc"    : 64, 
        "knl"    : 64
    }


    @staticmethod
    def checkPythonVersion():
        """check version. Python 3.3 required"""
        requiredVersion = (3,3)
        currentVersion  = sys.version_info
        if(requiredVersion > currentVersion):
            sys.exit("Requires Python 3.3 or newer. Abort.")
