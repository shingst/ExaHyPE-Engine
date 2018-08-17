
import sys

# add path to dependencies
from ..configuration import Configuration
sys.path.insert(1, Configuration.pathToCodegenerator)

import codegenerator

class CodegeneratorModel:

    # static list of dict, store all information on generated code (used for kernel registration)
    codegeneratorContextsList = []
    
    
    @staticmethod
    def generateCode(codegeneratorContext):
        # call the codegenerator with the given context
        codegeneratorController = codegenerator.Controller(codegeneratorContext)
        codegeneratorController.generateCode()
        # append the given context to the list
        CodegeneratorModel.codegeneratorContextsList.append(codegeneratorContext)