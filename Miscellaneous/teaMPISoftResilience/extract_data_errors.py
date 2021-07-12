#!/usr/bin/python

import sys
import csv
import re
#import matplotlib.pyplot as plt
import os
import math

# 220.799197   [i01r01c03s10.sng.lrz.de],rank:0, core:0, tid:0 info         exahype::runners::Runner::printTimeStepInfo(...)        step 11	team = 0	t_min          =0.009801
timestep_pattern=re.compile("( [0-9]*\.[0-9]*).*step.*team\ =\ ([0-9]).*t_min.*")
#162.815577   [i01r01c03s09.sng.lrz.de],rank:0, core:25, tid:0 info         exahype::offloading::ReplicationStatistics::printStatistics  team 1 spawned tasks = 156250 executed tasks = 68468 saved tasks =  78335 sent tasks = 68395 received tasks = 85139 received keys= 0 sent keys= 0 declined tasks = 0 late tasks = 3085
replication_pattern=re.compile("( [0-9]*\.[0-9]*).*ReplicationStatistics::printStatistics.*team\ ([0-9]).*saved tasks =  ([0-9]*) .*")

# 3.321646     [i01r09c04s10],rank:0, core:53, tid:1 error        exahype::reactive::ResilienceTools::overwriteDoubleIfActive() overwrite double value, pos = 906 old =-1.84735394850702130197364404377e-10 new = -1.00000000018473533813789799751 corresponds to relative error -230280961.399737566709518432617 corresponds to absolute error -1 max error indicator derivatives 1 max error indicator timestepsizes 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/reactive/ResilienceTools.cpp,line:197)
error_overwrite_pattern=re.compile("( [0-9]*\.[0-9]*).*ResilienceTools::overwrite.*IfActive.*old =(-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*new = (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*relative error (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*absolute error (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*max error indicator derivatives (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*).*max error indicator timestepsizes (-?[0-9]*\.?[0-9]*[e]?[-]?[0-9]*)")

# 3.322003     [i01r09c04s10],rank:0, core:53, tid:1 error        exahype::solvers::ADERDGSolver::computeErrorIndicator   Celldesc =15 errorIndicator derivative 290596004.635142624378204345703 errorIndicator time step size 0 errorIndicator admissibility 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver_MigratablePredictionJob.cpp,line:672)
error_indicator_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::computeErrorIndicator.*errorIndicator derivative ([0-9]*\.?[0-9]*) errorIndicator time step size ([0-9]*\.?[0-9]*)")

# 4.364058     [login03],rank:0, core:7, tid:0 error        exahype::solvers::ADERDGSolver::checkCellDescriptionAgainstOutcome soft error detected:  center[0] = 642.857143 center[1] = 357.142857 center[2] = 928.571429 time stamp = 0.034216 time step = 0.004277 element = 0 origin = 0 error indicator = 278.294152 isCorrupted = 0 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver.cpp,line:3885)
error_detect_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::checkCellDescriptionAgainstOutcome.*soft error detected")

# 4.364132     [login03],rank:0, core:7, tid:0 error        exahype::solvers::ADERDGSolver::correctWithOutcome()    Corrected an error in STP. (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver.cpp,line:3826)
error_correct_pattern=re.compile("( [0-9]*\.[0-9]*).*ADERDGSolver::correctWithOutcome.*Corrected an error in STP")

errors = [-1e-7,-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1e0,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]

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
  }
  print (filename)
  file = open(filename,"r")
  for line in file:
    #if "overwriteDouble" in line:
    print (line)
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

  #if(res["errors_corrected"]==1 and res["error"]>0):
  #  print (filename)
  if(res["error"]==0):
    print("Problem")

  print(res)
  return res


def printAccumResults(accumResults):
  print("max_error_indicator_derivative\tmax_error_indicator_ts\terror\truns\tdetected\tcorrected\tdetection_rate\tcorrection_rate")
  for res in accumResults:
    if res["runs"]>0:
      print(res["max_err_ind_der"],'\t', res["max_err_ind_ts"],'\t',res["error"],'\t',res["runs"],'\t',res["detected"],'\t',res["corrected"],'\t',res["detected"]/res["runs"],res["corrected"]/res["runs"])

 
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
                "runs" : float(res[4]) }
    results.append(res_dict)
  return results

def accumulateResults(results):
  unique_rel_errors = []
  [unique_rel_errors.append(res["error"]) for res in results if res["error"] not in unique_rel_errors]
  unique_rel_errors.sort()

  unique_error_indicators_der = []
  [unique_error_indicators_der.append(res["max_err_ind_der"]) for res in results if res["max_err_ind_der"] not in unique_error_indicators_der]
  unique_error_indicators_der.sort()
  print(unique_error_indicators_der)
  
  unique_error_indicators_ts = []
  [unique_error_indicators_ts.append(res["max_err_ind_ts"]) for res in results if res["max_err_ind_ts"] not in  unique_error_indicators_ts]
  unique_error_indicators_ts.sort()

  accumResults = []

  for max_error_indicator_der in unique_error_indicators_der:
    for max_error_indicator_ts in unique_error_indicators_ts:
      #for rel_error in unique_rel_errors:
      for error in errors:
        runs = 0
        detected = 0
        detected_corrected = 0

        #print(rel_error, max_error_indicator)

        for res in results:
          if(res["error"]==error and res["max_err_ind_der"]==max_error_indicator_der and res["max_err_ind_ts"]==max_error_indicator_ts):
            runs = runs+1
            if(res["errors_detected"]>0):
               detected=detected+1
            if(res["errors_detected"]==1 and res["errors_corrected"]==1):
               detected_corrected=detected_corrected+1

            if(runs==throw_away_results_at_run):
              break

        res = {"max_err_ind_der" : max_error_indicator_der,
               "max_err_ind_ts" : max_error_indicator_ts,
               "error" : error,
               "detected" : detected,
               "corrected" : detected_corrected,
               "runs" : runs
        }
        print (res)
        accumResults.append(res)
  return accumResults

if __name__=="__main__":
  dir_base = sys.argv[1]

  results = []

  for filename in os.listdir(dir_base):
    if filename.endswith(".err"):
      #print("reading", filename)
      result_tuple=parseErrorFile(dir_base+"/"+filename, 2)
      #print(result_tuple)
      results += [result_tuple]

  results = accumulateResults(results)
  printAccumResults(results)

  writeResults(results, "output.txt")
  results = readResults("output.txt")

  ##plt.savefig("error_sensitivity.pdf", bbox_inches='tight', dpi=300)
  #plt.show()
  

