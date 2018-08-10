import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required

    #TODO Add missing path
    # path from this file to ExaHyPE's root
    pathToExaHyPERoot          = os.path.join("..", "..", "..")
    
    # path to jinja2
    pathToJinja2               = os.path.join("..", "..", "dependencies", "jinja")
    
    # path to markupsafe
    pathToMarkupsafe           = os.path.join("..", "..", "dependencies", "markupsafe")
    
    alignmentPerArchitectures  = {
        "noarch" : 16,
        "snb"    : 32, 
        "hsw"    : 32, 
        "knc"    : 64, 
        "knl"    : 64
    }
