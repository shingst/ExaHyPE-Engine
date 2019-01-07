import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array

import re
import sys
import time

file = open(sys.argv[1], 'r')

ranks = int(sys.argv[2])
timestep_pattern = re.compile(".*step ([0-9]+).*t_min.*")
task_offload_pattern = re.compile(".*rank:([0-9]+).*resetRemainingTasksToOffload.*to rank ([0-9]+) ntasks ([0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")


cur_tasks_to_offload = np.zeros([ranks,ranks])
cur_blacklist_values = np.zeros([ranks,1])

#animation = sys.argv[3].lower() in ("yes", "true", "t", "1")
#if animation:
#  plt.ion()

current_step = -1

#fig, ax = plt.subplots()
#plt.show()

for line in file:
  m=timestep_pattern.match(line)
  if m:
    current_step = m.group(1)
    fig, ax = plt.subplots()
    ax.set_title("time step: "+current_step)
    print current_step
  
    for i in range(0,ranks):
      cur_tasks_to_offload[i][i]= -cur_blacklist_values[i]

    print cur_tasks_to_offload
    
    cur_tasks_to_offload_a = masked_array(cur_tasks_to_offload, cur_tasks_to_offload<0)
    cur_tasks_to_offload_b = masked_array(cur_tasks_to_offload, cur_tasks_to_offload>-0.5)
    
    pa = ax.matshow(cur_tasks_to_offload_a, cmap=plt.cm.Reds)
    cba = plt.colorbar(pa)
    pb = ax.matshow(cur_tasks_to_offload_b, cmap=plt.cm.gray, vmax=-0.5, vmin=-20)
    cbb = plt.colorbar(pb)

    for i in range(0,ranks):
      for j in range(0,ranks):
         val= cur_tasks_to_offload[j,i]
         ax.text(i, j, str(val), va='center', ha='center')

    plt.show()
   # if animation:
      #wait=raw_input("Press enter to continue")
   #   time.sleep(1)
      #plt.close()
  m=task_offload_pattern.match(line)
  if m:
    #print line
    cur_tasks_to_offload[int(m.group(1))][int(m.group(2))]=m.group(3)
    #print cur_tasks_to_offload
  m=blacklist_pattern.match(line)
  if m:
    #print line
    #print float(m.group(2))
    cur_blacklist_values[int(m.group(1))]=float(m.group(2))

file.close()
