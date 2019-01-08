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
task_offload_pattern = re.compile(".*rank:([0-9]+).*resetRemainingTasksToOffload.*to rank ([0-9]+) ntasks ([0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")


cur_tasks_to_offload = np.zeros([ranks,ranks])
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
    tmp_tasks_to_offload_a = masked_array(tmp_copy, tmp_copy<0)
    tmp_tasks_to_offload_b = masked_array(tmp_copy, tmp_copy>-0.5)
    
    pa = ax.matshow(tmp_tasks_to_offload_a, cmap=plt.cm.Reds)
    #cba = plt.colorbar(pa)
    pb = ax.matshow(tmp_tasks_to_offload_b, cmap=plt.cm.gray, vmax=-0.5, vmin=-20)
    #cbb = plt.colorbar(pb)

    if animate:
      ims_tmp.append(pa)
      ims_tmp.append(pb)

    for i in range(0,ranks):
      for j in range(0,ranks):
         val= cur_tasks_to_offload[j,i]
         txt = ax.text(i, j, str(val), va='center', ha='center')
         if animate:
           ims_tmp.append(txt)

    if animate:
      ims.append(ims_tmp) 
    if not animate:
      plt.show()
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


if animate:
  ani = animation.ArtistAnimation(fig, ims, interval=int(current_step), blit=True, repeat_delay=1000)
  ani.save("movie.mp4")
  #vid=ani.to_html5_video()
  #print vid
  plt.show()

file.close()
