#!/usr/bin/python

import sys
import csv
import re
#import matplotlib.pyplot as plt
import os
import math

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

#errors = [-1e-7,-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1e0,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
errors = [-1e-7,-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1e0,-1e1,-1e2,-1e3,-1e4,1e4,1e3,1e2,1e1,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]

throw_away_results_at_run = 100

def parseErrorFile(filename, teams):
  #dict={}
  #for i in range(teams):
  res={
      "old" : 0,
      "new" : 0,
      "error_indicator_der" : 0,
      "error_indicator_ts" : 0,
      "error" : 0,
      "errors_detected" : 0,
      "errors_corrected" : 0,
      "max_err_ind_der" : 0,
      "max_err_ind_ts" : 0,
      "run_number" : 0
  }
  print (filename)
  file = open(filename,"r")
  for line in file:
    #if "overwriteDouble" in line:
    #print (line)
    m = re.match(error_overwrite_pattern, line)
    if m:
     #print (line)
     res["old"] = float(m.group(2))
     res["new"] =float(m.group(3))
     rel_error = float(m.group(4))
     abs_error = float(m.group(5))
     #if(rel_error>0):
     #  error = rel_error
     #else:
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
  #if(res["errors_corrected"]==1 and res["error"]>0);
  #  print (filename)

  # additional check
  if(res["errors_corrected"] == 1):
    injected = 0
    healed = 0
    outfilename = filename[:-3]+"out"
    outfile = open(outfilename,"r")
    for line in outfile:
      m = re.match(stats_pattern, line)
      if m:
        injected = int(m.group(2))
        healed = int(m.group(3))
    if(injected!=healed):
      print("ERROR: was corrected on wrong team!")
      res["errors_corrected"] = 0
 
  if(res!=None and res["error"]==0):
    print("Problem: file ",filename, " does not contain valid output")
    os.remove(filename)
    return None

  #print(res)
  return res


def printAccumResults(accumResults):
  print("max_error_indicator_derivative\tmax_error_indicator_ts\terror\truns\tdetected\tcorrected\tdetection_rate\tcorrection_rate")#\trun_numbers")
  for res in accumResults:
    if res["runs"]>0:
      print(res["max_err_ind_der"],'\t',res["max_err_ind_ts"],'\t',res["error"],'\t',res["runs"],'\t',res["detected"],'\t',res["corrected"],'\t',res["detected"]/res["runs"],res["corrected"]/res["runs"])#,res["run_numbers"])

 
#def sortResults(result_tuples):
#  result_tuples.sort()
#  return result_tuples
#
def writeResults(accumResults, outputfile):
#  sortedResults = sortResults(result_tuples)
  f = open(outputfile,"w")
  f.write("max_error_indicator_derivative\tmax_error_indicator_ts\terror\truns\tdetected\tcorrected\tdetection_rate\tcorrection_rate\n")
  for res in accumResults:
    if res["runs"]>0:
      f.write(str(res["max_err_ind_der"])+'\t'+str(res["max_err_ind_ts"])+'\t'+str(res["error"])+'\t'+str(res["runs"])+'\t'+str(res["detected"])+'\t'+str(res["corrected"])+'\t'+str(res["detected"]/res["runs"])+'\t'+str(res["corrected"]/res["runs"])+'\n')

def readResults(resultsfile):
  f = open(resultsfile, "r")
  results = []
  # skip header
  f.readline() 
  for line in f:
    res = line.split('\t')
    res_dict  = {"max_err_ind_der" : float(res[0]),
                "max_err_ind_ts" : float(res[1]),
                "error" : float(res[2]),
                "detected" : float(res[4]),
                "corrected" : float(res[5]),
                "runs" : float(res[3]) }
    results.append(res_dict)
  return results

def accumulateResults(results):
  unique_rel_errors = []
  [unique_rel_errors.append(res["error"]) for res in results if res["error"] not in unique_rel_errors]
  unique_rel_errors.sort()

  unique_error_indicators_der = []
  [unique_error_indicators_der.append(res["max_err_ind_der"]) for res in results if res["max_err_ind_der"] not in unique_error_indicators_der]
  unique_error_indicators_der.sort()
  #print("unique error indicators", unique_error_indicators_der)
  
  unique_error_indicators_ts = []
  [unique_error_indicators_ts.append(res["max_err_ind_ts"]) for res in results if res["max_err_ind_ts"] not in  unique_error_indicators_ts]
  unique_error_indicators_ts.sort()
  #print("unique error indicators", unique_error_indicators_ts)

  accumResults = []

  for max_error_indicator_der in unique_error_indicators_der:
    for max_error_indicator_ts in unique_error_indicators_ts:
      #for rel_error in unique_rel_errors:
      for error in errors:
        runs = 0
        detected = 0
        detected_corrected = 0
        run_numbers = []

        for res in results:
          if(res["error"]==error and res["max_err_ind_der"]==max_error_indicator_der and res["max_err_ind_ts"]==max_error_indicator_ts):
            #print("accumulating ", res)
            runs = runs+1
            if(res["errors_detected"]>0):
               detected=detected+1
            if(res["errors_detected"]==1 and res["errors_corrected"]==1):
               detected_corrected=detected_corrected+1
 
            run_numbers += [int(res["run_number"])]

            if(runs==throw_away_results_at_run):
              break
        run_numbers.sort()
        res = {"max_err_ind_der" : max_error_indicator_der,
               "max_err_ind_ts" : max_error_indicator_ts,
               "error" : error,
               "detected" : detected,
               "corrected" : detected_corrected,
               "runs" : runs#,
    #           "run_numbers" : run_numbers
        }
        if(res["runs"]>0):
          print (res)
          accumResults.append(res)
  return accumResults

def computeAggregateAccuracy(accumResults):
  #unique_error_indicators_der = []
  #[unique_error_indicators_der.append(res["max_err_ind_der"]) for res in accumResults if res["max_err_ind_der"] not in unique_error_indicators_der]
  #unique_error_indicators_der.sort()
  #print("unique error indicators", unique_error_indicators_der)
  
  #unique_error_indicators_ts = []
  #[unique_error_indicators_ts.append(res["max_err_ind_ts"]) for res in accumResults if res["max_err_ind_ts"] not in  unique_error_indicators_ts]
  #unique_error_indicators_ts.sort()
  #print("unique error indicators", unique_error_indicators_ts)
 
  error_inds = []
  error_inds = [(res["max_err_ind_der"],res["max_err_ind_ts"]) for res in accumResults if (res["max_err_ind_der"],res["max_err_ind_ts"]) not in error_inds]
  error_inds = list(set(error_inds))

  accuracy_dict = {}
  for err_ind in error_inds:
    accuracy_dict[err_ind] = (0,0)

  for res in accumResults:
    key = (res["max_err_ind_der"],res["max_err_ind_ts"])
    old_tuple = accuracy_dict[key]
    print (res)
    new_tuple = (old_tuple[0]+1, old_tuple[1]+res["corrected"]/res["runs"])
    accuracy_dict[key] = new_tuple
  for key in accuracy_dict:
    cnt = accuracy_dict[key][0]
    sum_accuracy = accuracy_dict[key][1]
    accuracy_dict[key] = (cnt,sum_accuracy/cnt)
  
  return accuracy_dict


if __name__=="__main__":
  dir_base = sys.argv[1]

  results = []

  for filename in os.listdir(dir_base):
    if filename.endswith(".err"):
      #print("reading", filename)
      result_tuple=parseErrorFile(dir_base+"/"+filename, 2)
      #print(result_tuple)
      if(result_tuple!=None):
        results += [result_tuple]

  results = accumulateResults(results)
  #print(results)
  printAccumResults(results)

  writeResults(results, "output.txt")
  results = readResults("output.txt")

  accuracies = computeAggregateAccuracy(results)

  for key in sorted(accuracies):
    print (key, accuracies[key])
  ##plt.savefig("error_sensitivity.pdf", bbox_inches='tight', dpi=300)
  #plt.show()
  

