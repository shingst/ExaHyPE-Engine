import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required
    
    # absolute path to jinja2
    pathToJinja2      = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "jinja"))
    
    # absolute path to markupsafe (jinja2 dependency)
    pathToMarkupsafe  = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "markupsafe"))
    
    # absolute path to jsonschema
    pathToJSONSchema  = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "jsonschema"))
    
    # absolute path to attr (jsonschema dependency, inside attrs/src)
    pathToAttr        = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "attrs", "src"))
    
    # absolute path to pyrsistent (jsonschema dependency)
    pathToPyrsistent  = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "Submodules", "pyrsistent"))
    
    # absolute path to the schema file
    pathToSchema      = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "exahype-specfile.schema.json"))


    @staticmethod
    def checkPythonVersion():
        """check version. Python 3.3 required"""
        requiredVersion = (3,3)
        currentVersion  = sys.version_info
        if(requiredVersion > currentVersion):
            sys.exit("Requires Python 3.3 or newer. Abort.")
