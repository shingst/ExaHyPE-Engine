#!/usr/bin/python

import sys
import csv
import re
#import matplotlib.pyplot as plt
import os
import math

import distutils.util

# 220.799197   [i01r01c03s10.sng.lrz.de],rank:0, core:0, tid:0 info         exahype::runners::Runner::printTimeStepInfo(...)        step 11	team = 0	t_min          =0.009801
timestep_pattern=re.compile("( [0-9]*\.[0-9]*).*step.*team\ =\ ([0-9]).*t_min.*")

# 2.195630     [i04r11c03s05],rank:0, core:0, tid:0 info         exahype::reactive::JobTableStatistics::printStatistics   team 0 spawned tasks = 3773 executed tasks = 1909 double checked tasks = 0 soft errors injected = 1 detected soft errors = 0 limited tasks = 0 healed tasks = 0 saved tasks =  1864 sent tasks = 1909 received tasks = 1959 declined tasks = 0 late tasks = 72
stats_pattern=re.compile("( [0-9]*\.[0-9]*).*JobTableStatistics::printStatistics.*team\ 0.*soft errors injected = ([0-9]*) .*healed tasks = ([0-9]*)")

# 3.321646     [i01r09c04s10],rank:0, core:53, tid:1 error        exahype::reactive::ResilienceTools::overwriteDoubleIfActive() overwrite double value, pos = 906 old =-1.84735394850702130197364404377e-10 new = -1.00000000018473533813789799751 corresponds to relative error -230280961.399737566709518432617 corresponds to absolute error -1 max error indicator derivatives 1 max error indicator timestepsizes 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/reactive/ResilienceTools.cpp,line:197)
error_overwrite_pattern=re.compile("( [0-9]*\.[0-9]*).*ResilienceTools::overwrite.*IfActive.*old =(-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*new = (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*relative error (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*absolute error (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*max error indicator derivatives (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*max error indicator timestepsizes (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*)")

# 3.322003     [i01r09c04s10],rank:0, core:53, tid:1 error        exahype::solvers::ADERDGSolver::computeErrorIndicator   Celldesc =15 errorIndicator derivative 290596004.635142624378204345703 errorIndicator time step size 0 errorIndicator admissibility 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver_MigratablePredictionJob.cpp,line:672)
error_indicator_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::computeErrorIndicator.*errorIndicator derivative ([0-9]*\.?[0-9]*) errorIndicator time step size ([0-9]*\.?[0-9]*)")

# 4.364058     [login03],rank:0, core:7, tid:0 error        exahype::solvers::ADERDGSolver::checkCellDescriptionAgainstOutcome soft error detected:  center[0] = 642.857143 center[1] = 357.142857 center[2] = 928.571429 time stamp = 0.034216 time step = 0.004277 element = 0 origin = 0 error indicator = 278.294152 isCorrupted = 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver.cpp,line:3885)
error_detect_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::checkCellDescriptionAgainstOutcome.*soft error detected")

# 4.364132     [login03],rank:0, core:7, tid:0 error        exahype::solvers::ADERDGSolver::correctWithOutcome()    Corrected an error in STP. (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver.cpp,line:3826)
error_correct_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::correctWithOutcome.*Corrected an error in STP")

#sweep/parameters={"order": "5", "maximumMeshSize": "400", "maximumMeshDepth": "0", "kernels": "generic", "soft_error_generation": "migratable_stp_tasks_overwrite", "check_mechanism": "check_stps_with_low_confidence", "save_redundancy": "true", "task_sharing_mode": "task_sharing_resilience_correction", "check_time_steps": "true", "check_derivatives": "true", "check_admissibility": "true", "check_lazily": "true", "max_error_indicator_derivatives": "0", "max_error_indicator_timestepsizes": "0.03", "offloading": "none", "offloading_progress": "none", "offloading_CCP_temperature": "0.5", "offloading_diffusion_temperature": "0.5", "offloading_CCP_frequency": "0", "offloading_CCP_steps": "0", "offloading_update_temperature": "true", "offloading_increase_temp_threshold": "0", "abs_error": "-1", "rel_error": "0", "injection_rank": "0", "architecture": "skx", "dimension": "3", "bufferSize": "64", "timeSteps": "10", "timeStepping": "global"}
configuration_pattern=re.compile("sweep/parameters.*\"check_time_steps\": \"(true|false)\".*\"check_derivatives\": \"(true|false)\".*\"check_admissibility\": \"(true|false)\".*\"check_lazily\": \"(true|false)\".*")

#errors = [-1e-7,-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1e0,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
errors = [-1e-7,-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1e0,-1e1,-1e2,-1e3,-1e4,1e4,1e3,1e2,1e1,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]

throw_away_results_at_run = 100

def parseResult(filenamePrefix, teams):
  res={
      "old_value" : 0,
      "new_value" : 0,
      "check_adm" : False,
      "check_der" : False,
      "check_ts" : False,
      "check_lazy": False,
      "error_indicator_der" : 0,
      "error_indicator_ts" : 0,
      "error" : 0,
      "errors_detected" : 0,
      "errors_corrected" : 0,
      "max_err_ind_der" : 0,
      "max_err_ind_ts" : 0,
      "run_number" : 0
  }

  #parse error file
  errFilename = filenamePrefix+".err"
  #print(errFilename)
  err_file = open(errFilename,"r")
  for line in err_file:
    m = re.match(error_overwrite_pattern, line)
    if m:
     res["old"] = float(m.group(2))
     res["new"] =float(m.group(3))
     rel_error = float(m.group(4))
     abs_error = float(m.group(5))
     error = abs_error
     res["error"] = math.copysign(1,error) * 10 **math.ceil(math.log10(abs(error)))
     res["max_err_ind_der"] = float(m.group(6))
     res["max_err_ind_ts"] = float(m.group(7))
    m = re.match(error_indicator_pattern, line)
    if m:
     #print (line)
     res["error_indicator_der"] = float(m.group(2))
     res["error_indicator_ts"] = float(m.group(3))
    m = re.match(error_detect_pattern, line)
    if m:
     res["errors_detected"] +=1
    m = re.match(error_correct_pattern, line)
    if m:
     res["errors_corrected"] +=1

  res["run_number"] = int(filename.split("-")[-1][1:-4])

  # parse out file
  injected = 0
  healed = 0
  outFilename = filenamePrefix+".out"
  outfile = open(outFilename,"r")
  for line in outfile:
    m = re.match(stats_pattern, line)
    if m:
      injected = int(m.group(2))
      healed = int(m.group(3))
    m = re.match(configuration_pattern, line)
    if m:
      res["check_ts"] = (bool(distutils.util.strtobool(m.group(1))))
      res["check_der"] = (bool(distutils.util.strtobool(m.group(2))))
      res["check_adm"] = (bool(distutils.util.strtobool(m.group(3))))
      res["check_lazy"] = (bool(distutils.util.strtobool(m.group(4))))
  if(injected!=healed):
    print("ERROR: was not corrected appropriately:",outFilename)
    res["errors_corrected"] = 0
  elif (res["errors_corrected"]==1):
    print("Was corrected appropriately:",outFilename)

  if(res!=None and res["error"]==0):
    print("Problem: file ",errFilename, " does not contain valid output")
    os.remove(errFilename)
    return None

  return res

def accumResultToString(res):
  if res["runs"]>0:
    resStr = ""
    for k in res:
      resStr += str(res[k])+'\t'
    resStr += str(res["detected"]/res["runs"])+'\t'
    resStr += str(res["corrected"]/res["runs"])+'\t'
    return (resStr)
  else:
    return None

def printAccumResults(accumResults):
  print("check_adm\tcheck_ts\tcheck_der\tcheck_lazy\tmax_error_indicator_derivative\tmax_error_indicator_ts\terror\truns\tdetected\tcorrected\tdetection_rate\tcorrection_rate")#\trun_numbers")
  for res in accumResults:
    resstr = accumResultToString(res)
    if(resstr!=None):
      print(resstr)

def writeAccumResults(accumResults, outputfile):
#  sortedResults = sortResults(result_tuples)
  f = open(outputfile,"w")
  f.write("check_adm\tcheck_ts\tcheck_der\tcheck_lazy\tmax_error_indicator_derivative\tmax_error_indicator_ts\terror\truns\tdetected\tcorrected\tdetection_rate\tcorrection_rate\n")#\trun_numbers")
  for res in accumResults:
    resstr = accumResultToString(res)
    if(resstr!=None):
      f.write(resstr+"\n")
  f.close()

def readAccumResults(resultsfile):
  f = open(resultsfile, "r")
  results = []
  # skip header
  f.readline() 
  for line in f:
    res = line.split('\t')
    res_dict  = {
                "check_adm" : bool(distutils.util.strtobool(res[0])), 
                "check_ts" : bool(distutils.util.strtobool(res[1])), 
                "check_der" : bool(distutils.util.strtobool(res[2])), 
                "check_lazy" : bool(distutils.util.strtobool(res[3])), 
                "max_err_ind_der" : float(res[4]),
                "max_err_ind_ts" : float(res[5]),
                "error" : float(res[6]),
                "runs" : int(res[7]),
                "detected" : float(res[8]),
                "corrected" : float(res[9]),
                }
    results.append(res_dict)
  f.close()
  return results

def accumulateResults(results):
  #unique_rel_errors = []
  #[unique_rel_errors.append(res["error"]) for res in results if res["error"] not in unique_rel_errors]
  #unique_rel_errors.sort()

  #unique_error_indicators_der = []
  #[unique_error_indicators_der.append(res["max_err_ind_der"]) for res in results if res["max_err_ind_der"] not in unique_error_indicators_der]
  #unique_error_indicators_der.sort()
  #print("unique error indicators", unique_error_indicators_der)
  
  #unique_error_indicators_ts = []
  #[unique_error_indicators_ts.append(res["max_err_ind_ts"]) for res in results if res["max_err_ind_ts"] not in  unique_error_indicators_ts]
  #unique_error_indicators_ts.sort()
  #print("unique error indicators", unique_error_indicators_ts)

  unique_configs = [(res["max_err_ind_der"],res["max_err_ind_ts"], res["check_adm"], res["check_ts"], res["check_der"], res["check_lazy"]) for res in results]
  unique_configs = list(set(unique_configs))

  accumResults = []

  for config in unique_configs:
      for error in errors:
        runs = 0
        detected = 0
        detected_corrected = 0
        run_numbers = []

        for res in results:
          if(res["error"]==error and res["max_err_ind_der"]==config[0] and res["max_err_ind_ts"]==config[1] and
            res["check_adm"]==config[2] and res["check_ts"]==config[3] and res["check_der"]==config[4] and
            res["check_lazy"]==config[5]):
            runs = runs+1
            if(res["errors_detected"]>0):
               detected=detected+1
            if(res["errors_detected"]==1 and res["errors_corrected"]==1):
               detected_corrected=detected_corrected+1
 
            run_numbers += [int(res["run_number"])]

            if(runs==throw_away_results_at_run):
              break
        run_numbers.sort()
        res = {"check_adm" : config[2],
               "check_ts" : config[3],
               "check_der" : config[4],
               "check_lazy" : config[5],
               "max_err_ind_der" : config[0],
               "max_err_ind_ts" : config[1],
               "error" : error,
               "runs" : runs,
               "detected" : detected,
               "corrected" : detected_corrected,
    #           "run_numbers" : run_numbers
        }
        if(res["runs"]>0):
          #print (res)
          accumResults.append(res)
  return accumResults

def computeAvgSensitivity(accumResults):
  unique_configs = [(res["max_err_ind_der"],res["max_err_ind_ts"], res["check_adm"], res["check_ts"], res["check_der"],res["check_lazy"]) for res in accumResults]
  unique_configs = list(set(unique_configs))

  accuracy_res = []
  accuracy_dict = {}
  for config in unique_configs:
    accuracy_dict[config] = (0,0)

  for res in accumResults:
    key = (res["max_err_ind_der"],res["max_err_ind_ts"],res["check_adm"],res["check_ts"], res["check_der"], res["check_lazy"])
    old_tuple = accuracy_dict[key]
    #print (res)
    if(res["runs"]>0):
      new_tuple = (old_tuple[0]+1, old_tuple[1]+res["corrected"]/res["runs"])
    accuracy_dict[key] = new_tuple
  for key in accuracy_dict:
    cnt = accuracy_dict[key][0]
    if(cnt>0):
      sum_accuracy = accuracy_dict[key][1]
      accuracy_dict[key] = (cnt,sum_accuracy/cnt)
      res_element = {
			"check_adm" : key[2],
                     "check_ts" : key[3],
                     "check_der" : key[4],
                     "check_lazy" : key[5],
                     "max_err_ind_der" : key[0],
                     "max_err_ind_ts" : key[1],
                     "avg_sensitivity" : sum_accuracy/cnt,
                     "normalised_realtime_min" : 0,
                     "cnt" : cnt
                  }
      accuracy_res.append(res_element)
  return accuracy_res

def printSensitivityResults(accuracyResults):
  print("check_adm\tcheck_ts\tcheck_der\tcheck_lazy\tmax_err_ind_der\tmax_err_ind_ts\tavg_sensitivity\tnormalised_realtime_min\ttradeoff_ratio\tcnt")
  for res in accuracyResults:  
    if(res["normalised_realtime_min"]>0):
      factor = res["avg_sensitivity"]/res["normalised_realtime_min"]
    else:
      factor = 0
   
    print(res["check_adm" ],'\t',res["check_ts" ],\
         '\t',res["check_der"],'\t',res["check_lazy"],'\t',res["max_err_ind_der"],'\t',res["max_err_ind_ts"],'\t',res["avg_sensitivity"],'\t',res["normalised_realtime_min"],'\t',factor,'\t',res["cnt"])
    
if __name__=="__main__":
  dir_base = sys.argv[1]

  results = []

  for filename in os.listdir(dir_base):
    if filename.endswith(".err"):
      result_tuple=parseResult(dir_base+"/"+filename[0:-4], 2)
      #print(result_tuple)
      if(result_tuple!=None):
        results += [result_tuple]

  results = accumulateResults(results)
  printAccumResults(results)

  writeAccumResults(results, "output.txt")
  results = readAccumResults("output.txt")

  accuracies = computeAvgSensitivity(results)

  printSensitivityResults(accuracies) 

  

