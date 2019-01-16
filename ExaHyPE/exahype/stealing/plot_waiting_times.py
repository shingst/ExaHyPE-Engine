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
waiting_times_pattern = re.compile(".*rank:0.*updateLoadDistribution\(\) rank ([0-9]+) waiting for ([0-9]+) for rank ([0-9]+)")
#waiting_times_pattern = re.compile(".*rank:0.*updateLoadDistribution() rank ([0-9]+) waiting for ([0-9]+) for rank ([0-9]+)")

cur_waiting_times = np.zeros([ranks,ranks])

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

    pa = ax.matshow(cur_waiting_times, cmap=plt.cm.Reds)

    if animate:
      ims_tmp.append(pa)

    #for i in range(0,ranks):
    #  for j in range(0,ranks):
    #       blacklist_val = cur_tasks_to_offload[i,i]
    #       if blacklist_val<-0.5:
    #         print ('adding')
    #         txt = ax.text(i, j, str(blacklist_val), va='center', ha='center')
    #         if animate:
    #          ims_tmp.append(txt)

    #iif animate:
    #  ims.append(ims_tmp)
    if not animate:
      plt.show()
   
  m=waiting_times_pattern.match(line)
  if m:
    #print (line)
    cur_waiting_times[int(m.group(1))][int(m.group(3))]= int(m.group(2)) 

if animate:
  ani = animation.ArtistAnimation(fig, ims, interval=int(current_step), blit=True, repeat_delay=1000)
  #ani.save("movie.mp4")
  #vid=ani.to_html5_video()
  #print vid
  plt.show()

file.close()
