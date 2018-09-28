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
    """job_statistics.py:
    Plots a histogram from Peano's Job output

    run:

    ./job_statistics.py <file/folder> [output directory]

    file must contain output from an executable compiled with -DTBB_USE_THREADING_TOOLS
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
    stats = {}
    #  10.4573      info         no of background tasks[1]=19029
    ex = re.compile(r"background tasks\[(\d+)\]=(\d+)?")
    base = 10
    with open(fname, 'r') as f:
        for l in f.readlines():
            if "no of background tasks" in l:
                match = ex.search(l)
                if match:
                    num_tasks = int(match.group(1))
                    occurences = int(match.group(2))
                    if num_tasks == 0:
                        stats[0] = occurences
                    else:
                        stats_bin = int(math.log(num_tasks, base)) + 1
                        if stats_bin in stats:
                            stats[stats_bin] += occurences
                        else:
                            stats[stats_bin] = occurences
    if len(stats) == 0:
        print("WARNING: {} doesn't contain TBB_USE_THREADING_TOOLS output: no plot created".format(fname))
        return
    labels = []
    x = []
    y = []
    for k in range(min(stats.keys()), max(stats.keys())+1):
        x.append(k)
        if k in stats:
            y.append(stats[k])
        else:
            y.append(0)
        if k == 0:
            labels.append(0)
        else:
            labels.append(r"$<{}^{}$".format(base,k))
    plt.bar(x, y,tick_label=labels,log=True)
    plt.xlabel("Number of tasks executed by consumer")
    plt.ylabel("Frequency")

    if output_folder is None: output_folder="."
    out_file = os.path.join(output_folder, os.path.basename(fname)) + ".job-stats.pdf"
    plt.savefig(out_file, bbox_inches='tight')
    plt.clf()

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

    
