#!/usr/bin/python

import sys
import csv
import re
import matplotlib.pyplot as plt
import os

# 220.799197   [i01r01c03s10.sng.lrz.de],rank:0, core:0, tid:0 info         exahype::runners::Runner::printTimeStepInfo(...)        step 11	team = 0	t_min          =0.009801
timestep_pattern=re.compile("( [0-9]*\.[0-9]*).*step.*team\ =\ ([0-9]).*t_min.*")
#162.815577   [i01r01c03s09.sng.lrz.de],rank:0, core:25, tid:0 info         exahype::offloading::ReplicationStatistics::printStatistics  team 1 spawned tasks = 156250 executed tasks = 68468 saved tasks =  78335 sent tasks = 68395 received tasks = 85139 received keys= 0 sent keys= 0 declined tasks = 0 late tasks = 3085
replication_pattern=re.compile("( [0-9]*\.[0-9]*).*ReplicationStatistics::printStatistics.*team\ ([0-9]).*saved tasks =  ([0-9]*) .*")

# 0.536263     [login02],rank:0, core:2, tid:0 error        exahype::reactive::ResilienceTools::overwriteDoubleIfActive() overwrite double value, pos = 0 old =0 new = -0.001 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/reactive/ResilienceTools.cpp,line:131)
error_pattern=re.compile("( [0-9]*\.[0-9]*).*ResilienceTools::overwriteDoubleIfActive().*old =([0-9]*\.(0-9]*).*new = ([0-9]*\.[0-9]*) .*")
# 0.536293     [login02],rank:0, core:2, tid:0 error        exahype::solvers::ADERDGSolver::corruptIfActive         Has corrupted STP job  cellDescriptionIndex = 471 center[0] = 5.640000 center[1] = 0.600000 time stamp = 0.000000 time step = 0.077104 element = 0 origin = 0 confidence = 0.000000 isCorrupted = 1 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver_MigratablePredictionJob.cpp,line:379)
# 0.536309     [login02],rank:0, core:2, tid:0 error        exahype::solvers::ADERDGSolver::setConfidence           Celldesc =471 confidence 1 (file:/dss/dsshome1/02/di57zoh3/Codes/ExaHyPE-Engine/./ExaHyPE/exahype/solvers/ADERDGSolver_MigratablePredictionJob.cpp,line:624)

def parseFile(filename, teams):
  dict={}
  for i in range(teams):
    dict[i]={
      "timestamps" : [],
      "timestamps_delta" : [],
      "tasks_saved_accum" : [],
      "tasks_saved_delta" : []
    }

  file = open(filename,"r")
  for line in file:
    #if "ReplicationStatistics" in line:
    # print line
    m = re.match(timestep_pattern, line)
    if m:
     dict[int(m.group(2))]["timestamps"].append(float(m.group(1)))
    m = re.match(replication_pattern, line)
    if m:
     dict[int(m.group(2))]["tasks_saved_accum"].append(int(m.group(3)))

  for i in range(teams):
    for t in range(len(dict[i]["timestamps"])):
      if(i==0):
        otherteam=1
      else:
        otherteam=0
      dict[i]["timestamps_delta"].append(abs(dict[i]["timestamps"][t]-dict[otherteam]["timestamps"][t])) 

    previous_saved = 0 
    for t in range(len(dict[i]["tasks_saved_accum"])):
      current_saved=dict[i]["tasks_saved_accum"][t]
      dict[i]["tasks_saved_delta"].append(current_saved-previous_saved)
      previous_saved=current_saved
  return dict

def writeToFiles(outfile, teams, dict):
  fieldnames=['step','timestamps','timestamps_delta','tasks_saved_accum','tasks_saved_delta']
  for i in range(0,teams):
    with open(outfile+"_team_"+str(i), 'w') as file:
     writer=csv.DictWriter(file, fieldnames=fieldnames)
     writer.writeheader()
     for t in range(len(dict[i]['timestamps'])):
      if(len(dict[i]['tasks_saved_accum'])>0):
        writer.writerow({ 'step' : t, 'timestamps' : dict[i]['timestamps'][t], 'timestamps_delta' : dict[i]['timestamps_delta'][t], 'tasks_saved_accum' : dict[i]['tasks_saved_accum'][t],
                      'tasks_saved_delta' : dict[i]['tasks_saved_delta'][t] })
      else:
        writer.writerow({ 'step' : t, 'timestamps' : dict[i]['timestamps'][t], 'timestamps_delta' : dict[i]['timestamps_delta'][t], 'tasks_saved_accum' : 0,
                      'tasks_saved_delta' : 0 })

def plotDivergence(x, data, title, marker, colour, facecolors):
  handle=plt.scatter(x, data, label=title, marker=marker, color=colour, s=32,facecolors=facecolors)
  return handle

def plotSavedTasks(x, data, title, marker, colour):
  handle=plt.scatter(x, data, label=title, marker=marker, facecolors='none',edgecolors=colour, s=16, alpha=0.5)
  return handle


if __name__=="__main__":
  dir_base = sys.argv[1]
  dir_sharing = sys.argv[2]

  dicts_no_sharing = []
  dicts_sharing = []

  for run in range(1,21):
    for filename in os.listdir(dir_base):
      if filename.endswith("-r"+str(run)+".out"):
        dict_no_sharing=parseFile(dir_base+"/"+filename, 2)
        writeToFiles("no_sharing_r"+str(run), 2, dict_no_sharing)
        dicts_no_sharing.append(dict_no_sharing)

    for filename in os.listdir(dir_sharing):
      if filename.endswith("-r"+str(run)+".out"):
        dict_sharing=parseFile(dir_sharing+"/"+filename, 2)
        writeToFiles("sharing_r"+str(run), 2, dict_sharing)
        dicts_sharing.append(dict_sharing)

    #print run  
    #handle_no_sharing = plotDivergence(range(0,len(dicts_no_sharing[run-1][1]['timestamps_delta'])),dicts_no_sharing[run-1][1]['timestamps_delta'], 'no task sharing', 's','g','none')
    #handle_sharing = plotDivergence(range(0,len(dicts_sharing[run-1][1]['timestamps_delta'])),dicts_sharing[run-1][1]['timestamps_delta'], 'task sharing', 'o','b','b')
  
  #plt.xlim(left=-1)
  #plt.ylim(bottom=-15, top=70)
  #plt.xticks([0,5,10,15,20])
  #plt.xlabel('time step')
  #plt.ylabel('divergence between team 0 and team 1')
  #plt.legend((handle_no_sharing, handle_sharing),('task sharing disabled','task sharing enabled'), loc='center right')
  #plt.tight_layout()
  #plt.savefig("divergence.pdf", dpi=300)
  #plt.show()

  for run in range(1,21):
    handle_t1 = plotSavedTasks(range(0,len(dicts_sharing[run-1][1]['tasks_saved_delta'])), dicts_sharing[run-1][1]['tasks_saved_delta'], 'team 1 (disturbed)', 'o', 'k')
    handle_t0 = plotSavedTasks(range(0,len(dicts_sharing[run-1][0]['tasks_saved_delta'])), dicts_sharing[run-1][0]['tasks_saved_delta'], 'team 0', 's','r')
 
  avgs_t1 = []
  avgs_t0 = []

  for l in range(0,len(dicts_sharing[0][1]['tasks_saved_delta'])):   
    sum_t1 = 0
    sum_t0 = 0
    for run in range(1,21):
      sum_t1 += dicts_sharing[run-1][1]['tasks_saved_delta'][l]
      sum_t0 += dicts_sharing[run-1][0]['tasks_saved_delta'][l]
    avgs_t1.append(sum_t1/20)
    avgs_t0.append(sum_t0/20) 
 
  handle_avg_1 = plt.plot(avgs_t1, color='k',marker='o',linewidth=3, markersize=10)
  handle_avg_0 = plt.plot(avgs_t0, color='r',marker='s',linewidth=3, markersize=10)

  plt.xlim(left=-1)
  plt.xlim(right=21)
  plt.xticks([0,5,10,15,20])
  plt.ylim(bottom=0)
  plt.xlabel('time step')
  plt.ylabel('number of reused tasks per team per time step')
  plt.legend((handle_t1, handle_t0, handle_avg_1[0], handle_avg_0[0]),('team 1 (disturbed)', 'team 0','team 1 (disturbed) - average', 'team 0 - average'), loc='upper left')
  plt.tight_layout()
  plt.savefig("saved_tasks.pdf", bbox_inches='tight', dpi=300)
  plt.show()
  

