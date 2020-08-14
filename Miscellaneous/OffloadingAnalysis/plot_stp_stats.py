import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array

import re
import sys
import time

file = open(sys.argv[1], 'r')

ranks = int(sys.argv[2])
timesteps = int(sys.argv[3])
timestep_pattern = re.compile(".([0-9]+\.[0-9]+).*step ([0-9]+).*t_min.*")
blacklist_pattern = re.compile(".*blacklist value for rank ([0-9]+):([0-9]+\.[0-9]+)")
# 226.719163   [lns08.hpc.itc.rwth-aachen.de],rank:3, core:17, tid:0 info         exahype::stealing::AggressiveHybridDistributor::printOffloadingStatistics() time per STP  0.00216728
stp_pattern = re.compile(".*rank:([0-9]+).*printOffloadingStatistics.*time per STP  ([0-9]+\.[0-9]+)")


current_step = -1
last_timestamp = 0
duration = -1
duration_arr = []

x=[[0.0 for _ in range(timesteps)] for _ in range(ranks)]
cnt=[0 for _ in range(ranks)]

for line in file:
  m=stp_pattern.match(line)
  if m:
    print (line)
    rank  = int(m.group(1))
    time = float(m.group(2))
    if(cnt[rank]>=timesteps): 
       continue

    x[rank][cnt[rank]]=time
    cnt[rank]=cnt[rank]+1

print x


ts=range(0,timesteps)

fig,ax1 = plt.subplots()
for i in range(1,ranks):
 plt.scatter(ts,x[i])
plt.show()
#ax1.plot(duration_arr, 'b-x')
#ax1.set_xlabel('time step')
#ax1.set_ylabel('time step duration (wall clock)', color='b')
#ax1.set_ylim(bottom=0)

#ax2 = ax1.twinx()
#ax2.plot(tasksoffloaded_arr, 'r-x')
#ax2.set_ylabel('number of tasks offloaded', color='r')
#
#fig.tight_layout()
#plt.savefig("timestepstats.pdf")
file.close()

