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
# Generates solver header and source files
#


from .abstractModelBaseClass import AbstractModelBaseClass
from .codegeneratorModel import CodegeneratorModel


class SolverModel(AbstractModelBaseClass):
    
    def generateADERDGSolverFiles(self):
        solverTemplates = {
            "user"    : [ (self.context["solver"]+".h"   , "solvers/MinimalADERDGSolverHeader.template"),
                          (self.context["solver"]+".cpp" , "solvers/EmptyADERDGSolverImplementation.template") ],
            "generic" : [ (self.context["solver"]+".h"   , "solvers/ADERDGSolverHeader.template"),
                          (self.context["solver"]+".cpp" , "solvers/ADERDGSolverInCUserCode.template") ],
        }
        solverTemplates["optimised"] = solverTemplates["generic"]
        
        abstractSolverTemplates  = { 
            "user"      :  [],
            "generic"   :  [ (self.context["abstractSolver"]+".cpp" , "solvers/AbstractGenericADERDGSolverImplementation.template"),
                             (self.context["abstractSolver"]+".h"   , "solvers/AbstractGenericADERDGSolverHeader.template") ],
            "optimised" :  [ (self.context["abstractSolver"]+".cpp" , "solvers/AbstractOptimisedADERDGSolverImplementation.template"),
                             (self.context["abstractSolver"]+".h"   , "solvers/AbstractOptimisedADERDGSolverHeader.template") ],
        }
        
        implementation = self.context["implementation"]
        
        result = []
        for filePath,template in solverTemplates.get(implementation,[]):
          result.append(self.render(template,filePath,overwrite=False))
        for filePath,template in abstractSolverTemplates.get(implementation,[]):
          result.append(self.render(template,filePath))
        
        if implementation == "optimised":
            CodegeneratorModel.generateCode(self.context)
        
        return filter(lambda x: x is not None, result) # return generated files as list, None from not overwrite is filtered out

    def generateFiniteVolumesSolverFiles(self):
        solverTemplates = {
            "user"    : [ (self.context["solver"]+".h"   , "solvers/MinimalFiniteVolumesSolverHeader.template"),
                          (self.context["solver"]+".cpp" , "solvers/EmptyFiniteVolumesSolverImplementation.template") ],
            "generic" : [ (self.context["solver"]+".h"   , "solvers/FiniteVolumesHeader.template"),
                          (self.context["solver"]+".cpp" , "solvers/FiniteVolumesInCUserCode.template") ],
        }
        
        abstractSolverTemplates  = { 
            "generic"   :  [ (self.context["abstractSolver"]+".cpp" , "solvers/AbstractGenericFiniteVolumesSolverImplementation.template"),
                             (self.context["abstractSolver"]+".h"   , "solvers/AbstractGenericFiniteVolumesSolverHeader.template") ]
        }
       
        implementation = self.context["implementation"]
        
        if implementation=="optimised": # TODO
            print("ERROR: optimised FV kernels not available yet.",file=sys.stderr)
            raise
            
        result = []
        for filePath,template in solverTemplates.get(implementation,[]):
          result.append(self.render(template,filePath,overwrite=False))
        for filePath,template in abstractSolverTemplates.get(implementation,[]):
          result.append(self.render(template,filePath))
        
        return filter(lambda x: x is not None, result) # return generated files as list, None from not overwrite is filtered out

    def generateLimitingADERDGSolverFiles(self):
        solverTemplates = {  }
        abstractSolverTemplates  = { 
            "generic"   :  [ (self.context["abstractSolver"]+".cpp" , "solvers/AbstractGenericLimiterSolverImplementation.template"),
                             (self.context["abstractSolver"]+".h"   , "solvers/AbstractGenericLimiterSolverHeader.template") ],
            "optimised" :  [ (self.context["abstractSolver"]+".cpp" , "solvers/AbstractOptimisedLimiterSolverImplementation.template"),
                             (self.context["abstractSolver"]+".h"   , "solvers/AbstractOptimisedLimiterSolverHeader.template") ],
        }
        
        implementation = self.context["implementation"]
        
        if implementation=="user": # TODO
            print("ERROR: optimised FV kernels not available yet.",file=sys.stderr)
            raise
        
        result = []
        for filePath,template in abstractSolverTemplates.get(implementation,[]):
          result.append(self.render(template,filePath))
        
        return result # return generated files as list
        
    def generateCode(self):
        generators = { 
          "ADER-DG"          : self.generateADERDGSolverFiles,
          "Finite-Volumes"   : self.generateFiniteVolumesSolverFiles,
          "Limiting-ADER-DG" : self.generateLimitingADERDGSolverFiles
        }
       
        return generators[self.context["solverType"]]()
