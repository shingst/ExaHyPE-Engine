import os
import jinja2

import helper

import generator
import plotter

##
# Generates the solver files
class SolverGenerator(generator.Generator):
    
    def __init__(self,spec,specFile,verbose):
        generator.Generator.__init__(self,spec,specFile,verbose)
        self.plotterGenerator = plotter.PlotterGenerator(spec,specFile,verbose)
    
    def writeSolverFiles(self,solverMap,abstractSolverMap,implementation,context):
        for filePath in solverMap.get(implementation,{}):
            generator.Generator.render_template(self,solverMap[implementation][filePath],context,filePath,"Generated user solver file",False)
        for filePath in abstractSolverMap.get(implementation,{}):
            generator.Generator.render_template(self,abstractSolverMap[implementation][filePath],context,filePath,"Generated abstract solver file",True)
    
    def generateADERDGSolverFiles(self,context):
        solverMap = {
            "user"    :  { context["solver"]+".h"               : "solvers/MinimalADERDGSolverHeader.template", 
                           context["solver"]+".cpp"             : "solvers/EmptyADERDGSolverImplementation.template" },
            "generic" : {  context["solver"]+".h"               : "solvers/ADERDGSolverHeader.template", 
                           context["solver"]+".cpp"             : "solvers/ADERDGSolverInCUserCode.template"},
        }
        solverMap["optimised"] = solverMap["generic"]
        
        abstractSolverMap  = { 
            "user"      :  { },
            "generic"   :  { context["abstractSolver"]+".cpp"    : "solvers/AbstractGenericADERDGSolverImplementation.template", 
                             context["abstractSolver"]+".h"      : "solvers/AbstractGenericADERDGSolverHeader.template" },
            "optimised" :  { context["abstractSolver"]+".cpp"    : "solvers/AbstractOptimisedADERDGSolverImplementation.template", 
                             context["abstractSolver"]+".h"      : "solvers/AbstractOptimisedADERDGSolverHeader.template" }
        }
        
        implementation = context["implementation"]
        try:
            self.writeSolverFiles(solverMap,abstractSolverMap,implementation,context)
        except Exception:
            raise
        if implementation=="optimised":
            print("ERROR: not implemented yet '"+filePath+"'",file=sys.stderr)
            sys.exit(-1)

    ##
    # Generate user solver and abstract solver header and source files for an FV solver.
    #
    # @note Does not overwrite user solver files if they already exist.
    #
    def generateFVSolverFiles(self,context):
        solverMap = {
            "user"    : { context["solver"]+".h"   : "solvers/MinimalFiniteVolumesSolverHeader.template", 
                          context["solver"]+".cpp" : "solvers/EmptyFiniteVolumesSolverImplementation.template" },
            "generic" : { context["solver"]+".h"   : "solvers/FiniteVolumesHeader.template", 
                          context["solver"]+".cpp" : "solvers/FiniteVolumesInCUserCode.template"},
        }
        
        abstractSolverMap  = { 
            "generic"   :  { context["abstractSolver"]+".cpp" : "solvers/AbstractGenericFiniteVolumesSolverImplementation.template", 
                             context["abstractSolver"]+".h"   : "solvers/AbstractGenericFiniteVolumesSolverHeader.template" }
        }
        
        implementation = context["implementation"]
        try:
            self.writeSolverFiles(solverMap,abstractSolverMap,implementation,context)
        except Exception:
            raise
        
        if implementation=="optimised":
            print("ERROR: not implemented yet '"+filePath+"'",file=sys.stderr)
            sys.exit(-1)
    
    def generateLimitingADERDGSolverFiles(self,context):
        solverMap = {  }
        abstractSolverMap  = { 
            "generic"   :  { context["abstractSolver"]+".cpp" : "solvers/AbstractGenericLimiterSolverImplementation.template", 
                             context["abstractSolver"]+".h"   : "solvers/AbstractGenericLimiterSolverHeader.template" },
            "optimised" :  { context["abstractSolver"]+".cpp" : "solvers/AbstractOptimisedLimiterSolverImplementation.template", 
                             context["abstractSolver"]+".h"   : "solvers/AbstractOptimisedLimiterSolverHeader.template" }
        }
        implementation = solver["fv_kernel"].get("implementation","generic")
        try:
            self.writeSolverFiles(solverMap,abstractSolverMap,implementation,context)
        except Exception:
            raise
