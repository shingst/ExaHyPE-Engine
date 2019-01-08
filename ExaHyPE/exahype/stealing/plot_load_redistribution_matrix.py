import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.animation as animation

import re
import sys
import time

file = open(sys.argv[1], 'r')

ranks = int(sys.argv[2])
timestep_pattern = re.compile(".*step ([0-9]+).*t_min.*")
task_offload_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* target tasks to rank ([0-9]+) ntasks ([0-9]+) not offloaded ([0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")


cur_tasks_to_offload = np.zeros([ranks,ranks])
cur_tasks_not_offloaded = np.zeros([ranks,ranks])
cur_blacklist_values = np.zeros([ranks,1])

animate = sys.argv[3].lower() in ("yes", "true", "t", "1")
print (animate)
#if animate:
#  plt.ion()

current_step = -1

if animate:
  fig, ax = plt.subplots()
  ims = []


for line in file:
  m=timestep_pattern.match(line)
  if m:
    current_step = m.group(1)
    if not animate:
      fig, ax = plt.subplots()
      ax.set_title("time step: "+current_step)
    if animate:
      ims_tmp = []
    print (current_step)
  
    for i in range(0,ranks):
      cur_tasks_to_offload[i][i]= -cur_blacklist_values[i]

    print (cur_tasks_to_offload)
   
    tmp_copy = cur_tasks_to_offload.copy()

    for i in range(0,ranks):
      for j in range(0,ranks):
         if(i!=j):
           tasks_to_offload = cur_tasks_to_offload[j,i]
           tasks_not_offloaded = cur_tasks_not_offloaded[j,i]
           tasks_offloaded = tasks_to_offload - tasks_not_offloaded
           tmp_copy[j][i] = tasks_offloaded

    tmp_tasks_to_offload_a = masked_array(tmp_copy, tmp_copy<0)
    tmp_tasks_to_offload_b = masked_array(tmp_copy, tmp_copy>-0.5)
    
    pa = ax.matshow(tmp_tasks_to_offload_a, cmap=plt.cm.Reds, vmin=0, vmax=20000)
    #cba = plt.colorbar(pa)
    pb = ax.matshow(tmp_tasks_to_offload_b, cmap=plt.cm.gray, vmax=-0.5, vmin=-20)
    #cbb = plt.colorbar(pb)

    if animate:
      ims_tmp.append(pa)
      ims_tmp.append(pb)
    
    for i in range(0,ranks):
      for j in range(0,ranks):
         if(i!=j):
           tasks_to_offload = cur_tasks_to_offload[j,i]
           tasks_not_offloaded = cur_tasks_not_offloaded[j,i]
           tasks_offloaded = tasks_to_offload - tasks_not_offloaded
           txt = ax.text(i, j, str(tasks_offloaded)+"/"+str(tasks_to_offload), va='center', ha='center')
           if animate:
             ims_tmp.append(txt)
         else:
           blacklist_val = cur_tasks_to_offload[i,i]
           if blacklist_val<-0.5:
             print ('adding')
             txt = ax.text(i, j, str(blacklist_val), va='center', ha='center')
             if animate:
               ims_tmp.append(txt)

    if animate:
      ims.append(ims_tmp) 
    if not animate:
      plt.show()
  m=task_offload_pattern.match(line)
  if m:
    #print line
    cur_tasks_to_offload[int(m.group(1))][int(m.group(2))] = m.group(3)
    cur_tasks_not_offloaded[int(m.group(1))][int(m.group(2))]= m.group(4)
    #print cur_tasks_to_offload
  m=blacklist_pattern.match(line)
  if m:
    print line
    #print float(m.group(2))
    cur_blacklist_values[int(m.group(1))]=float(m.group(2))


if animate:
  ani = animation.ArtistAnimation(fig, ims, interval=int(current_step), blit=True, repeat_delay=1000)
  #ani.save("movie.mp4")
  #vid=ani.to_html5_video()
  #print vid
  plt.show()

file.close()
