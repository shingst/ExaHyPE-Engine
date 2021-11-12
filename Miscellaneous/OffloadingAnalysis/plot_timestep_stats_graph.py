import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array
from graphviz import Digraph

import re
import sys
import time

file = open(sys.argv[1], 'r')

ranks = int(sys.argv[2])
timestep_pattern = re.compile(".([0-9]+\.[0-9]+).*step ([0-9]+).*t_min.*")
task_offload_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.* target tasks to rank ([0-9]+) ntasks ([0-9]+) not offloaded ([0-9]+).*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")


dot = Digraph()

for i in range(0,ranks):
  dot.node(str(i), str(i))


current_step = -1
last_timestamp = 0
duration = -1
duration_arr = []


for line in file:
  m=timestep_pattern.match(line)
  if m:
    #print line
    current_step = m.group(2)
    if current_step>0:
      duration = float(m.group(1))-last_timestamp
    last_timestamp = float(m.group(1))
    duration_arr.append(duration)
  
  
    dot.render('test_graph'+str(current_step))    

    dot = Digraph()

    for i in range(0,ranks):
      dot.node(str(i), str(i))
  m=task_offload_pattern.match(line)
  if m:
    print (line)
    tasksoffloaded = int(m.group(3)) - int(m.group(4))
    if tasksoffloaded>0:
      dot.edge(m.group(1),m.group(2), label=str(tasksoffloaded)+"/"+m.group(3))


file.close()
