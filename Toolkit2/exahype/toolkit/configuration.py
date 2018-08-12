import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required

    #TODO Add missing path
    # absolute path to ExaHyPE's root (we need absolute paths in generated Makefile)
    pathToExaHyPERoot          = os.path.abspath(os.path.join(__file__,"..", "..", "..", ".."))
    
    # absolute path to jinja2
    pathToJinja2               = os.path.abspath(os.path.join(__file__,"..", "..", "dependencies", "jinja"))
    
    # absolute path to markupsafe
    pathToMarkupsafe           = os.path.abspath(os.path.join(__file__,"..", "..", "dependencies", "markupsafe"))
    
    alignmentPerArchitectures  = {
        "noarch" : 16,
        "snb"    : 32, 
        "hsw"    : 32, 
        "knc"    : 64, 
        "knl"    : 64
    }
