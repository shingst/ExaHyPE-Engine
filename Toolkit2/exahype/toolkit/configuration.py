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
    pathToExaHyPERoot   = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    
    # absolute path to jinja2
    pathToJinja2        = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "jinja"))
    
    # absolute path to markupsafe (jinja2 dependency)
    pathToMarkupsafe    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "markupsafe"))
    
    # absolute path to the codegenerator module
    pathToCodegenerator = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "CodeGenerator"))
    
    # absolute path to the specfile module
    pathToSpecfiles    = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
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
