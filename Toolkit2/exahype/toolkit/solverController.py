import os
import copy

from .models import *
from .toolkitHelper import ToolkitHelper

class SolverController:

    def __init__(self, solverSpec, baseContext):
        self.solverSpec = solverSpec
        self.baseContext = baseContext


    def run(self, logger):
        for i,solver in enumerate(self.solverSpec):
            logger.info("Generating solver[%d] = %s..." % (i, solver["name"]))
            solverFiles = []
            if solver["type"]=="ADER-DG":
                model       = solverModel.SolverModel(self.buildADERDGSolverContext(solver))
                solverFiles = model.generateCode()
            elif solver["type"]=="Finite-Volumes":
                model       = solverModel.SolverModel(self.buildFVSolverContext(solver))
                solverFiles = model.generateCode()
            elif solver["type"]=="Limiting-ADER-DG":
                aderdgContext = self.buildADERDGSolverContext(solver)
                fvContext     = self.buildFVSolverContext(solver)
                context       = self.buildLimitingADERDGSolverContext(solver)
                # modifications
                fvContext["solver"]         = context["FVSolver"]
                fvContext["solverType"]     = "Finite-Volumes"
                fvContext["abstractSolver"] = context["FVAbstractSolver"]
                fvContext["patchSize"]      = 2 * aderdgContext["order"] + 1
                
                aderdgContext["solver"]                 = context["ADERDGSolver"]
                aderdgContext["solverType"]             = "ADER-DG"
                aderdgContext["abstractSolver"]         = context["ADERDGAbstractSolver"]
                aderdgContext["numberOfDMPObservables"] = context["numberOfDMPObservables"]
                aderdgContext["ghostLayerWidth"]        = fvContext["ghostLayerWidth"]
                
                # generate all solver files
                model        = solverModel.SolverModel(context)
                solverFiles  = model.generateCode()
                model        = solverModel.SolverModel(fvContext)
                solverFiles += model.generateCode()
                model        = solverModel.SolverModel(aderdgContext)
                solverFiles += model.generateCode()
            
            for path in solverFiles:
                logger.info("Generated '"+path+"'")
            for j,plotter in enumerate(solver.get("plotters",[])):
                logger.info("Generating plotter[%d] = %s for solver[%d] = %s" % (j, plotter["name"], i, solver["name"]))
                model = plotterModel.PlotterModel(self.buildPlotterContext(solver,plotter))
                for path in model.generateCode():
                    logger.info("Generated '"+path+"'")
                        


    def buildBaseSolverContext(self, solver):
        context = copy.copy(self.baseContext)
        
        context["solverType"]     = solver["type"]
        context["solver"]         = solver["name"]
        context["abstractSolver"] = "Abstract"+context["solver"]
        
        nVar          = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"variables"))
        nParam        = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"material_parameters"))
        nGlobalObs    = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"global_observables"))
        nPointSources = solver["point_sources"] if type(solver.get("point_sources",[])) is int else len(solver.get("point_sources",[]))
        
        context["numberOfVariables"]          = nVar
        context["numberOfMaterialParameters"] = nParam
        context["numberOfGlobalObservables"]  = nGlobalObs
        context["numberOfPointSources"]       = nPointSources
        
        # variables access class
        context["variablesMap"]  = ToolkitHelper.parse_variables(solver,"variables")
        if nParam>0:
            context["variablesMap"] += ToolkitHelper.parse_variables(solver,"material_parameters")
        context["variablesMapSize"] = len(context["variablesMap"])
        
        context["range_0_nVar"]          = range(0,nVar)
        context["range_0_nVarParam"]     = range(0,nVar+nParam)
        context["range_0_nGlobalObs"]    = range(0,nGlobalObs)    # nGlobalObs might be 0
        context["range_0_nPointSources"] = range(0,nPointSources) # nPointSources might be 0
        
        context["namingSchemes"]={} # TODO read from spec
        return context


    def buildADERDGSolverContext(self, solver):
        context = self.buildBaseSolverContext(solver)
        context.update(self.buildADERDGKernelContext(solver["aderdg_kernel"]))
        context.update(self.buildKernelOptimizationContext(solver["aderdg_kernel"], context))
        
        context["order"]                  = solver["order"]
        context["numberOfDMPObservables"] = -1 # overwrite if called from LimitingADERDGSolver creation
        
        return context
        
    def buildLimitingADERDGSolverContext(self, solver):
        context = self.buildBaseSolverContext(solver)
        
        context["order"]                  = solver["order"]
        context["numberOfDMPObservables"] = solver["limiter"]["dmp_observables"]
        context["implementation"]         = solver["limiter"].get("implementation","generic")
        
        context["ADERDGSolver"]         = solver["name"]+"_ADERDG"
        context["FVSolver"]             = solver["name"]+"_FV"
        context["ADERDGAbstractSolver"] = "Abstract"+solver["name"]+"_ADERDG"
        context["FVAbstractSolver"]     = "Abstract"+solver["name"]+"_FV"
        
        return context
    
        
    def buildFVSolverContext(self, solver):
        context = self.buildBaseSolverContext(solver)
        context.update(self.buildFVKernelContext(solver["fv_kernel"]))
        
        context["patchSize"] = solver.get("patch_size",-1) # overwrite if called from LimitingADERDGSolver creation
        
        return context
    

    
    def buildADERDGKernelContext(self, kernel):
        context = {}
        context["implementation"]          = kernel.get("implementation","generic")
        context["useMaxPicardIterations"]  = "true" if kernel.get("space_time_predictor",{}).get("maxpicarditer",0)!=0 else "false"
        context["maxPicardIterations"]     = kernel.get("space_time_predictor",{}).get("maxpicarditer",0)
        context["tempVarsOnStack"]         = kernel.get("allocate_temporary_arrays","heap")=="stack" 
        context["patchwiseAdjust"]         = kernel.get("adjust_solution","pointwise")=="patchwise" 
        context["language"]                = kernel.get("language","C").lower()
        context["basis"]                   = kernel.get("basis","Legendre").lower()
        context["isLinear"]                = not kernel.get("nonlinear",True)
        context["isNonlinear"]             = kernel.get("nonlinear",True)
        context["linearOrNonlinear"]       = "Linear" if context["isLinear"] else "Nonlinear"
        context["isFortran"]               = kernel.get("language",False)=="Fortran" 
        context["useCERK"]                 = kernel.get("space_time_predictor",{}).get("cerkguess",False)
        context["noTimeAveraging"]         = "true" if kernel.get("space_time_predictor",{}).get("notimeavg",False) else "false"
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context
        
    def buildFVKernelContext(self,kernel):
        context = {}
        ghostLayerWidth = { "godunov" : 1, "musclhancock" : 2 }
        context["ghostLayerWidth"]   = ghostLayerWidth[kernel["scheme"]]
        context["finiteVolumesType"] = kernel["scheme"]
        context["implementation"]    = kernel.get("implementation","generic")
        context["tempVarsOnStack"]   = kernel.get("allocate_temporary_arrays","heap")=="stack" 
        context["patchwiseAdjust"]   = kernel.get("adjust_solution","pointwise")=="patchwise" 
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context
    
    def buildKernelTermsContext(self,terms):
        context = {}
        for term in ["flux","source","ncp","point_sources","material_parameters"]:
            option = term.replace("_s","S").replace("_p","P").replace("ncp","NCP")
            option = "use%s%s" % ( option[0].upper(), option[1:] )
            context[option] = term in terms
            context["%s_s" % option] = "true" if context[option] else "false"
        return context

        
    def buildPlotterContext(self,solver,plotter):
        context = self.buildBaseSolverContext(solver)

        context["plotter"]         = plotter["name"]
        context["writtenUnknowns"] = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"variables"))
        context["plotterType"]     = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
        
        context["headerPath"] = os.path.join(context["plotterSubDirectory"],(context["plotter"]+".h"))
        if context["plotterSubDirectory"] != "":
            context["outputPath"] = os.path.join(context["outputPath"], context["plotterSubDirectory"])
            
        return context
        
    def buildKernelOptimizationContext(self, kernel, solverContext):
        optimizations = kernel.get("optimised_kernel_debugging",[]) + kernel.get("optimised_terms",[])
        context = {}
        context["optKernelPath"]      = os.path.join("kernels", solverContext["project"] + "_" + solverContext["solver"])
        context["optNamespace"]       = solverContext["project"] + "::" + solverContext["solver"] + "_kernels::aderdg"
        context["useConverter"]       = "converter"       in optimizations
        context["countFlops"]         = "flops"           in optimizations
        context["useFluxVect"]        = "fluxvect"        in optimizations
        context["useFusedSource"]     = "fusedsource"     in optimizations
        context["useFusedSourceVect"] = "fusedsourcevect" in optimizations
        context["useSourceVect"]      = "fusedsourcevect" in optimizations
        context["useNCPVect"]         = "fusedsourcevect" in optimizations

        return context
        
