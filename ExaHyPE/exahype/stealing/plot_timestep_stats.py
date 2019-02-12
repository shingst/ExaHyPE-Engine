import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array

import re
import sys
import time

file = open(sys.argv[1], 'r')

ranks = int(sys.argv[2])
timestep_pattern = re.compile(".([0-9]+\.[0-9]+).*step ([0-9]+).*t_min.*")
task_offload_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* target tasks to rank ([0-9]+) ntasks ([0-9]+) not offloaded ([0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")

#cur_tasks_to_offload = np.zeros([ranks,ranks])
#cur_tasks_not_offloaded = np.zeros([ranks,ranks])
#cur_blacklist_values = np.zeros([ranks,1])

#tasks_to_offload = np.zeros([ranks,ranks])
#tasks_to_offload2 = np.zeros([ranks,ranks])
#lines_read_for_rank = np.zeros([ranks,1])
#lines_read_for_rank2 = np.zeros([ranks,1])

#lines_read_critical = 0
#tasksoffloaded_critical = 0

current_step = -1
last_timestamp = 0
duration = -1
duration_arr = []

tasksoffloaded = 0
tasksoffloaded_arr = []

for line in file:
  m=timestep_pattern.match(line)
  if m:
    print (line)
    current_step = int(m.group(2))
    if current_step>0:
      duration = float(m.group(1))-last_timestamp
    last_timestamp = float(m.group(1))
    if current_step>0:
      print (duration)
      duration_arr.append(duration)
  
    tasksoffloaded_arr.append(tasksoffloaded)
    tasksoffloaded = 0
  m=task_offload_pattern.match(line)
  if m:
    print (line)
    tasksoffloaded += int(m.group(3)) - int(m.group(4))
  #  if(tasks>0):
  #    tasksoffloaded_critical += tasks
  #    lines_read_critical += 1
  #    if(lines_read_critical==ranks-1):
  #      lines_read_critical = 0
  #      tasksoffloaded_arr.append(tasksoffloaded_critical)
  #      tasksoffloaded_critical = 0
  #  if(lines_read_for_rank[int(m.group(1))]>=ranks-1):
  #    tasks_to_offload2[int(m.group(1))][int(m.group(2))] = int(m.group(3))-int(m.group(4))
  #    lines_read_for_rank2[int(m.group(1))] = lines_read_for_rank2[int(m.group(1))]+1
  #  else:
  #    tasks_to_offload[int(m.group(1))][int(m.group(2))] = int(m.group(3))-int(m.group(4))
  #    lines_read_for_rank[int(m.group(1))] = lines_read_for_rank[int(m.group(1))]+1
  #    print (lines_read_for_rank == ranks-1).all()
  #    if (lines_read_for_rank == ranks-1).all():
  #       lines_read_for_rank = lines_read_for_rank2.copy()
  #       tasksoffloaded_arr.append(np.sum(tasks_to_offload))
  #       tasks_to_offload = tasks_to_offload2.copy()
  #       tasks_to_offload2 = np.zeros([ranks, 1])
  #       lines_read_for_rank2 = np.zeros([ranks, 1]) 
#
 #     print ("buffer1",lines_read_for_rank)
#      print ("buffer2",lines_read_for_rank2)
  #m=blacklist_pattern.match(line)
  #if m:
  #  #print line
  #  #print float(m.group(2))
  #  cur_blacklist_values[int(m.group(1))]=float(m.group(2))

print (duration_arr)
print (tasksoffloaded_arr)

fig,ax1 = plt.subplots()

ax1.plot(duration_arr, 'b-x')
ax1.set_xlabel('time step')
ax1.set_ylabel('time step duration (wall clock)', color='b')
ax1.set_ylim(bottom=0)

ax2 = ax1.twinx()
ax2.plot(tasksoffloaded_arr, 'r-x')
ax2.set_ylabel('number of tasks offloaded', color='r')

fig.tight_layout()
plt.savefig("timestepstats.pdf")
file.close()
