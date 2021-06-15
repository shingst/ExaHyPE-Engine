#!/usr/bin/python

import sys
import csv
import re
import matplotlib.pyplot as plt
import os

# 220.799197   [i01r01c03s10.sng.lrz.de],rank:0, core:0, tid:0 info         exahype::runners::Runner::printTimeStepInfo(...)        step 11	team = 0	t_min          =0.009801
timestep_pattern=re.compile("( [0-9]*\.[0-9]*).*step.*team\ =\ ([0-9]).*t_min.*([0-9]*\.[0-9]*)")
#162.815577   [i01r01c03s09.sng.lrz.de],rank:0, core:25, tid:0 info         exahype::offloading::ReplicationStatistics::printStatistics  team 1 spawned tasks = 156250 executed tasks = 68468 saved tasks =  78335 sent tasks = 68395 received tasks = 85139 received keys= 0 sent keys= 0 declined tasks = 0 late tasks = 3085
# 11.185617    [philipp-ThinkPad-Edge-E540],rank:0, core:2, tid:0 info         exahype::reactive::JobTableStatistics::printStatistics   team 1 spawned tasks = 16250 executed tasks = 15657 double checked tasks = 14971 soft errors injected = 0 detected soft errors = 688 limited tasks = 0 healed tasks = 0 saved tasks =  0 sent tasks = 15579 received tasks = 15599 declined tasks = 0 late tasks = 0

detected_errors_pattern=re.compile("( [0-9]*\.[0-9]*).*JobTableStatistics::printStatistics.*team ([0-9]).*soft errors injected = ([0-9]*).*detected soft errors = ([0-9]*).*limited tasks = ([0-9]*).*")


def parseFile(filename, teams):
  dict={}
  for i in range(teams):
    dict[i]={
      "timesteps" :  [],
      "timestamps" : [],
      "timestamps_delta" : [],
      "detected_errors_accum" : [],
      "detected_errors_delta" : [],
      "injected_errors_accum" : [],
      "injected_errors_delta" : [],
      "limited_tasks_accum" : [],
      "limited_tasks_delta": []
    }

  file = open(filename,"r")
  for line in file:
    #if "ReplicationStatistics" in line:
    #print (line)
    m = re.match(timestep_pattern, line)
    if m:
     dict[int(m.group(2))]["timestamps"].append(float(m.group(1)))
     dict[int(m.group(2))]["timesteps"].append(float(m.group(3)))
    m = re.match(detected_errors_pattern, line)
    if m:
     dict[int(m.group(2))]["detected_errors_accum"].append(int(m.group(4)))
     dict[int(m.group(2))]["injected_errors_accum"].append(int(m.group(3)))
     dict[int(m.group(2))]["limited_tasks_accum"].append(int(m.group(5)))

  for i in range(teams):
    for t in range(len(dict[i]["timestamps"])):
      if(i==0):
        otherteam=1
      else:
        otherteam=0
      dict[i]["timestamps_delta"].append(abs(dict[i]["timestamps"][t]-dict[otherteam]["timestamps"][t])) 

    previous_detected = 0 
    for t in range(len(dict[i]["detected_errors_accum"])):
      current_detected=dict[i]["detected_errors_accum"][t]
      dict[i]["detected_errors_delta"].append(current_detected-previous_detected)
      previous_detected=current_detected
    
    previous_injected = 0 
    for t in range(len(dict[i]["injected_errors_accum"])):
      current_injected=dict[i]["injected_errors_accum"][t]
      dict[i]["injected_errors_delta"].append(current_injected-previous_injected)
      previous_injected=current_injected

    previous_limited = 0
    for t in range(len(dict[i]["limited_tasks_accum"])):
      current_limited=dict[i]["limited_tasks_accum"][t]
      dict[i]["limited_tasks_delta"].append(current_limited-previous_limited)
      previous_limited=current_limited
  return dict

def writeToFiles(outfile, teams, dict):
  fieldnames=['step','timestamps','timestamps_delta','detected_errors_accum','detected_errors_delta']
  for i in range(0,teams):
    with open(outfile+"_team_"+str(i), 'w') as file:
     writer=csv.DictWriter(file, fieldnames=fieldnames)
     writer.writeheader()
     for t in range(len(dict[i]['timestamps'])):
      if(len(dict[i]['detected_errors_accum'])>0):
        writer.writerow({ 'step' : t, 'timestamps' : dict[i]['timestamps'][t], 'timestamps_delta' : dict[i]['timestamps_delta'][t], 'detected_errors_accum' : dict[i]['detected_errors_accum'][t],
                      'detected_errors_delta' : dict[i]['detected_errors_delta'][t] })
      else:
        writer.writerow({ 'step' : t, 'timestamps' : dict[i]['timestamps'][t], 'timestamps_delta' : dict[i]['timestamps_delta'][t], 'detected_errors_accum' : 0,
                      'detected_errors_delta' : 0 })

def plotDetectedErrors(ax, x, data, label, marker, colour):
  handle=ax.plot(x, data, label=label, marker=marker, markersize=7, markerfacecolor='None', color=colour)
  return handle

def plotTimeStepSizesDelta(ax, x, data, label, marker, colour):
  handle=ax.plot(x, data, label=label, marker=marker, color=colour,alpha=0.5)
  return handle


if __name__=="__main__":
  inputFilename = sys.argv[1]
  result = parseFile(inputFilename, 2)

  fig, ax = plt.subplots()
  #print (result)
  #handle_detected = plotDetectedErrors(ax, range(0,len(result[0]['timestamps'])), result[0]['detected_errors_accum'], "with injected error - accumulated detected errors ", "s", "r")

  #ax2 = ax.twinx()
  #timestamps_div = [result[1]['timesteps'][i]-result[0]['timesteps'][i] for i in range(0,len(result[0]['timesteps']))]
  #print (timestamps_div)
  #handle_divergence = plotTimeStepSizesDelta(ax2,range(0, len(result[1]['timesteps'])), timestamps_div, "with injected error - divergence of time step sizes", "o", "b")

  handle_limited = plotTimeStepSizesDelta(ax,range(0, len(result[0]['limited_tasks_accum'])), result[0]['limited_tasks_accum'], "with injected error - limited tasks", "x", "r")
 
  for i in range(len(result[0]['injected_errors_delta'])):
    if (result[0]['injected_errors_delta'][i]>=1):
       ax.axvline(x=i, color='black', linestyle='--')
       print("error at ",i)

  for i in range(len(result[0]['limited_tasks_delta'])):
    if (result[0]['limited_tasks_delta'][i]>=1):
       ax.axvline(x=i, color='red', linestyle='--')
       print("limiter at ",i)
  

  #plt.xlim(left=-1)
  #plt.xlim(right=21)
  #plt.xticks([0,5,10,15,20])
  #plt.ylim(bottom=0)
  ax.set_xlabel('Time step')
  #ax.set_ylabel('Cumulative Number of Detected Errors on Team 0')
  ax.set_ylabel('Cumulative Number of Limited Tasks on Team 0')
  ax.set_ylim(bottom=0)
 
  #ax2.set_ylabel('Divergence of Timestep Sizes between Teams')
  #ax2.set_ylim([0,1e-05])
  #ax2.set_ylabel('Limited Tasks')
  #ax2.set_ylim([0,1e-05])

  # added these three lines
  lns = handle_limited
  #lns = handle_detected+handle_divergence
  labs = [l.get_label() for l in lns]
  ax.legend(lns, labs, loc=0)
  #ax.legend()
  #ax2.legend()
  #plt.legend((handle_t1, handle_t0, handle_avg_1[0], handle_avg_0[0]),('team 1 (disturbed)', 'team 0','team 1 (disturbed) - average', 'team 0 - average'), loc='upper left')
  plt.tight_layout()
  plt.savefig("limited_tasks_small_error.pdf", bbox_inches='tight', dpi=300)
  plt.show()
  

