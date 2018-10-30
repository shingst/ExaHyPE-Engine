#!/usr/bin/env python3

import sys,os
import matplotlib
import matplotlib.pyplot as plt
import math
import re

def haveToPrintHelpMessage(argv):
    """
    Check if we have to print a help message.
    """
    result = parseArgument(argv,1)==None
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def printHelpMessage():
    info = \
    """stealing_statistics.py:
    Plots various statistics from ExaHyPE's stealing profiler output

    run:

    ./stealing_statistics.py <file/folder> <ranks> [output directory]

    NOTE: check the output is not filtered out by a log filer
    """
    print(info) # correctly indented
    sys.exit()

def parseArgument(argv,i):
    if i<len(argv):
        return argv[i]
    else:
        return None

def plotStatistics(fname, output_folder):
    """
    79.827520    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler Stealing statistics for rank 1:
    79.827557    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   spawned tasks: 91125
    79.827578    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   performance updates: 161961
    79.827596    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   stealing decisions: 16412
    79.827610    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   target offloaded tasks:
    79.827626    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 0 : 20059
    79.827640    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 1 : 0
    79.827653    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 2 : 19472
    79.827666    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 3 : 19475
    79.827680    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total target offloaded tasks: 59006
    79.827694    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   offloaded tasks:
    79.827710    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 0 : 14250
    79.827738    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 1 : 0
    79.827753    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 2 : 12636
    79.827768    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     to rank: 3 : 12636
    79.827783    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total offloaded tasks: 39522
    79.827804    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   offloadable tasks that failed threshold requirement: 22901
    79.827819    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   received tasks:
    79.827840    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     from rank: 0 : 0
    79.827856    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     from rank: 1 : 0
    79.827870    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     from rank: 2 : 0
    79.827885    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler     from rank: 3 : 0
    79.827899    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total received tasks: 0
    79.827918    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total computation time: 251.630178
    79.827940    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total handling time: 0.000000
    79.827959    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total useful communication time: 38.718052
    79.828141    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total idle communication time: 6.643303
    79.828157    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total wait time for enclave tasks: 0.000000
    79.828180    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total wait time for skeleton tasks: 0.000002
    79.828195    [i25r03c01s10],rank:1, core:7, tid:24 info         exahype::stealing::StealingProfiler::exahype::stealing::StealingProfiler   total wait time for remote tasks: 0.000000
    """
  
    offloadedTasks = {}
    spawnedTasks = {}
    receivedTasks = {}
    totalCompTime = {}
    totalWaitTime = {}

    basename = os.path.basename(fname)
    if(output_folder is None): output_folder="."

    with open(fname, 'r') as f:
        for l in f.readlines():
            if "spawned tasks" in l:
                ex = re.compile(r"rank:(\d).*spawned tasks: (\d+)?")
                match = ex.search(l)
                if match:
                    rank = int(match.group(1))
                    num_tasks = int(match.group(2))
                    spawnedTasks[rank]= num_tasks
            if "total received tasks" in l:
                ex = re.compile(r"rank:(\d).*total received tasks: (\d+)?")
                match = ex.search(l)
                if match:
                    rank = int(match.group(1))
                    received = int(match.group(2)) 
                    receivedTasks[rank] = received
            if "total offloaded tasks" in l:
                ex = re.compile(r"rank:(\d).*total offloaded tasks: (\d+)?")
                match = ex.search(l)
                if match:
                    rank = int(match.group(1))
                    offloaded = int(match.group(2)) 
                    offloadedTasks[rank] = offloaded
            if "total computation time" in l:
                ex = re.compile(r"rank:(\d).*total computation time: (\d+\.?\d*)")
                match = ex.search(l)
                if match:
                    rank = int(match.group(1))
                    time = float(match.group(2))
                    totalCompTime[rank] = time
            if "total wait time" in l:
                ex = re.compile(r"rank:(\d).*total wait time.*: (\d+\.?\d*)")
                match = ex.search(l)
                if match:
                    rank = int(match.group(1))
                    time = float(match.group(2))
                    if rank in totalWaitTime:
                      totalWaitTime[rank]+=time
                    else: 
                      totalWaitTime[rank]=time

  
    BarWidth= 0.5
    labels = []
    x = []
    y = []
  
    for k in range(min(totalWaitTime.keys()), max(totalWaitTime.keys())+1):
        x.append(k)
        y.append(totalWaitTime[k])
        labels.append(str(k))
  
    plt.bar(x,y, label="total wait time", tick_label=labels, width=BarWidth)
    plt.xlabel("rank")
    plt.ylabel("time [s]")
    plt.title("total wait time")
    out_file= output_folder+"/"+basename+"_wait.pdf"
    #print (out_file)
    plt.savefig(out_file, bbox_inches='tight')
    #plt.show()
    plt.close()
 
    BarWidth= 0.5
    labels = []
    x = []
    y = []
  
    for k in range(min(offloadedTasks.keys()), max(offloadedTasks.keys())+1):
        x.append(k)
        y.append(offloadedTasks[k])
        labels.append(str(k))
  
    plt.bar(x,y, label="offloaded tasks", tick_label=labels, width=BarWidth)
    plt.xlabel("rank")
    plt.ylabel("tasks")
    plt.title("total offloaded tasks")
    out_file= output_folder+"/"+basename+"_offloaded.pdf"
    #print (out_file)
    plt.savefig(out_file, bbox_inches='tight')
    #plt.show()
    plt.close()

    BarWidth= 0.5
    labels = []
    x = []
    y = []
    
    for k in range(min(receivedTasks.keys()), max(receivedTasks.keys())+1):
        x.append(k)
        y.append(receivedTasks[k])
        labels.append(str(k))
  
    plt.bar(x,y, label="received tasks", tick_label=labels, width=BarWidth)
    plt.xlabel("rank")
    plt.ylabel("tasks")
    plt.title("total received tasks")
    out_file= output_folder+"/"+basename+"_received.pdf"
    plt.savefig(out_file, bbox_inches='tight')
    #plt.show()
    plt.close()

    BarWidth= 0.5
    labels = []
    x = []
    y = []
  
    for k in range(min(spawnedTasks.keys()), max(spawnedTasks.keys())+1):
        x.append(k)
        y.append(spawnedTasks[k])
        labels.append(str(k))
  
    plt.bar(x,y, label="spawned tasks", tick_label=labels, width=BarWidth)
    plt.xlabel("rank")
    plt.ylabel("tasks")
    plt.title("total spawned tasks")
    out_file= output_folder+"/"+basename+"_spawned.pdf"
    plt.savefig(out_file, bbox_inches='tight')
    #plt.show()
    plt.close()

    BarWidth= 0.5
    labels = []
    x = []
    y = []
    
    for k in range(min(totalCompTime.keys()), max(totalCompTime.keys())+1):
        x.append(k)
        y.append(totalCompTime[k])
        labels.append(str(k))
  
    plt.bar(x,y, label="total computation time", tick_label=labels, width=BarWidth)
    plt.xlabel("rank")
    plt.ylabel("time [s]")
    out_file= output_folder+"/"+basename+"_computation_time.pdf"
    plt.title ("total computation time")
    plt.savefig(out_file, bbox_inches='tight')
    #plt.show()
    plt.close()

if __name__ == "__main__":
    if haveToPrintHelpMessage(sys.argv):
        printHelpMessage()
    
    path_name = parseArgument(sys.argv,1)
    output_folder = parseArgument(sys.argv,2)

    if os.path.isfile(path_name):
        plotStatistics(path_name,output_folder)
    elif os.path.isdir(path_name):
        for f in os.listdir(path_name):
            plotStatistics(os.path.join(path_name,f), output_folder)
    else:
        print("ERROR: {} not recognized as a file or folder".format(sys.argv[1]))

    
