import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array
from pygraphviz import *
import matplotlib.image as mpimg
import matplotlib.animation as animation

import re
import sys,getopt
import time

timestep_pattern = re.compile(".([0-9]+\.[0-9]+).*step ([0-9]+).*t_min.*")
task_offload_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* target tasks to rank ([0-9]+) ntasks ([0-9]+) not offloaded ([0-9]+).*")
temperatureCCP_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* temperature value CCP ([0-9]+\.?[0-9]*).*")
temperatureDiffusion_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* temperature value diffusion ([0-9]+\.[0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")
waiting_times_pattern = re.compile(".*rank:0.*printWaitingTimes\(\) rank ([0-9]+) waiting for ([0-9]+\.[0-9]+) for rank ([0-9]+)")
critical_rank_pattern = re.compile(".*rank:([0-9]+).*updateLoadDistribution.*\(\).*optimal victim: ([0-9]+) critical rank:([0-9]+)")
critical_rank_pattern2 = re.compile(".*rank:([0-9]+).*updateLoadDistribution\(\).*current critical rank:([0-9]+)")

filename = ''
ranks = -1
animate = 0 

plot_waiting_times = 14
plot_tasks = 1
plot_ccp_temp = 1
plot_diffusion_temp = 1
plot_bval = 1

try: 
  opts, args = getopt.getopt(sys.argv[1:],"hf:r:a",["nowaittimes","notasks","noccptemp","nodifftemp","nobval"])
except getopt.GetoptError:
  print ('plot_waiting_times_graph.py -f input -r nranks [-a]')
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print ('plot_waiting_times_graph.py -f input -r nranks [-a]')
    sys.exit(2)
  elif opt in ("-f"):
    print (arg)
    filename = arg
  elif opt in ("-r"):
    ranks = int(arg)
  elif opt in ("-a"):
    animate = 1
  elif opt in ("--nowaittimes"):
    plot_waiting_times = 0
  elif opt in ("--notasks"):
    plot_tasks = 0
  elif opt in ("--noccptemp"):
    plot_ccp_temp = 0
  elif opt in ("--nodifftemp"):
    plot_diffusion_temp = 0
  elif opt in ("--nobval"):
    plot_bval = 0

file = open(filename, 'r')
#animate = sys.argv[3].lower() in ("yes", "true", "t", "1")

if animate:
  fig, ax = plt.subplots()
  ims =[]

dot = AGraph(strict=False,directed=True, style="filled",label="timestep: 0")

for i in range(0,ranks):
  dot.add_node(str(i))

zero_threshold = 0.020000

current_step = -1
last_timestamp = 0
duration = -1
duration_arr = []
current_waiting_times = np.zeros([ranks,ranks])

for line in file:
  m=timestep_pattern.match(line)
  if m:
    #print line
    current_step = int(m.group(2))
    if current_step>0:
      duration = float(m.group(1))-last_timestamp
    last_timestamp = float(m.group(1))
    duration_arr.append(duration)
  
    print ("step: ", current_step)
    output_file = "test_graph%03d.png" % current_step
    dot.draw(output_file, prog="dot")    

    if animate:
      im = plt.imshow(mpimg.imread(output_file))
      #plt.show()
      ims.append([im])

    dot = AGraph(strict=False,directed=True,style="filled", label="timestep: "+str(current_step+1))

    for i in range(0,ranks):
      dot.add_node(str(i), label=str(i))
  m=task_offload_pattern.match(line)
  if m:
    #print (line)
    tasksoffloaded = int(m.group(3)) - int(m.group(4))
    if tasksoffloaded>0 and plot_tasks:
      dot.add_edge(m.group(1),m.group(2), label=str(tasksoffloaded)+"/"+m.group(3))
  m=waiting_times_pattern.match(line)
  if m:
    src = int(m.group(1))
    dest = int(m.group(3))
    time = float(m.group(2))
    if time>zero_threshold and plot_waiting_times:
     current_waiting_times[int(m.group(1))][int(m.group(3))]= float(m.group(2))
     dot.add_edge(src, dest, label=str(time),fontcolor="red", color="red")
  m = critical_rank_pattern.match(line)
  if m:
    #print (line)
    printing_rank = int(m.group(1)) 
    optimal_victim = int(m.group(2))
    critical_rank = int(m.group(3))
    if(printing_rank==critical_rank):
      n=dot.get_node(critical_rank)
      n.attr['style']='filled'
      n.attr['fillcolor']='red'
      opt= dot.get_node(optimal_victim)
      opt.attr['style'] = 'filled'
      opt.attr['fillcolor']='green'
  m = critical_rank_pattern2.match(line)
  if m:
    print (line)
    printing_rank = int(m.group(1)) 
    critical_rank = int(m.group(2))
    if(printing_rank==critical_rank):
      n=dot.get_node(critical_rank)
      n.attr['style']='filled'
      n.attr['fillcolor']='red'
  m=blacklist_pattern.match(line)
  if m:
    #print line
    #print float(m.group(2))
    if float(m.group(2))>0.5:
      n = dot.get_node(int(m.group(1)))
      if plot_bval:
        n.attr['label']= m.group(1)+" bval="+m.group(2)
      n.attr['fillcolor'] = 'grey'
      n.attr['style'] = 'filled'
    #cur_blacklist_values[int(m.group(1))]=float(m.group(2))
  m=temperatureCCP_pattern.match(line)
  if m and current_step>0 and plot_ccp_temp:
    n = dot.get_node(int(m.group(1)))
    if(n.attr['label'] ==None):
      n.attr['label']=m.group(1)
    print (str(n.attr['label'])+" tempCCP="+m.group(2))
    n.attr['label']= str(n.attr['label'])+ " tempCCP="+m.group(2)
  m=temperatureDiffusion_pattern.match(line)
  if m and current_step>0 and plot_diffusion_temp:
    n = dot.get_node(int(m.group(1)))
    if(n.attr['label'] ==None):
      n.attr['label']=m.group(1)
    print (str(n.attr['label'])+" tempDif="+m.group(2))
    n.attr['label']= str(n.attr['label'])+ " tempDif="+m.group(2)

if animate:
    ani = animation.ArtistAnimation(fig, ims, interval=5000, blit=True, repeat_delay=1000)
    #ani.save("movie.avi")
    plt.show()

file.close()
