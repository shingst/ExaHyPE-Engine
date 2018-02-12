#!/usr/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# FROM: https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True}
plt.rcParams.update(params) 

import operator
import csv

'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Compare fused vs. nonfused time stepping.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Compare fused vs. nonfused time stepping.
'''

def column(matrix, i):
    return [row[i] for row in matrix]

def plotTimeStepComparison(table,identifier,nodes,tasks,cores,sharedMem,dimension,cells):
    '''
    Compare fused vs. nonfused time stepping via a bar chart.
    
    assume row layout:
        0        1    2     3       4        5      6          7               8           9       10           11               12        
    Identifier|Order|CC|Kernels|Algorithm|Adapter|Nodes|Tasks (per Node)|Cores (per Task)|Mode|Iterations|User Time (Total)|CPU Time (Total)
    '''
    adapters = {}
    adapters['nonfused'] = [ 
                                'BroadcastGlobalDataAndMergeNeighbourMessages',
                                'SolutionUpdate',
                                'Prediction'
                           ]
    adapters['fused']    = [ 
                                'FusedTimeStep'
                           ]
    
    labelsNonfused        = [ 
                       'Riemann',
                       'Correction',
                       'Prediction'
                    ]
    hatchNonfused = [ "////", "----",  "\\\\\\\\" ]

    
    # Plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    datafile = open(table, 'r')
    next(datafile) # skip header
    data = list(csv.reader(datafile,delimiter='&'))
    datafile.close()
    
    orders  = sorted(list(set(column(data,1))),key=str) 
    
    centers = range(0,len(orders))
    
    N    = 2
    yMin = 10.0**20
    yMax = 0.0
 
    foundAny = False
    for i,order in enumerate(orders):
        numberOfQuadraturePoints = ( cells*(int(order)+1)**dimension ) # TIME PER Q
        
        for j,algorithm in enumerate(['nonfused', 'fused']):
            # assume row layout:
            #      0        1    2     3       4        5      6          7               8           9       10           11               12        
            # Identifier|Order|CC|Kernels|Algorithm|Adapter|Nodes|Tasks (per Node)|Cores (per Task)|Mode|Iterations|User Time (Total)|CPU Time (Total)
            filtered = list(filter(lambda x: 
                  x[0]==identifier and x[1]==order and 
                  x[4]==algorithm and 
                  x[5] in adapters[algorithm] and
                  x[6]==nodes and x[7]==tasks and x[8]==cores and
                  x[9]==sharedMem,
                  data))
            
            totalRunTime = 0
            handles = []

            if len(filtered) and algorithm=='nonfused':
                foundAny=True
                
                for k,adapter in enumerate(adapters[algorithm]):
                    row           = list(filter(lambda x: x[5]==adapter,filtered))[0]
                    iterations    = float (row[10])
                    runtime       = float (row[11]) / numberOfQuadraturePoints / iterations
                    if i==0:
                        ax.bar(centers[i]-0.2,runtime,width=0.4,color='1.0',bottom=totalRunTime,align='center',log=False,hatch=hatchNonfused[k],label=labelsNonfused[k])
                    else:
                        ax.bar(centers[i]-0.2,runtime,width=0.4,color='1.0',bottom=totalRunTime,align='center',log=False,hatch=hatchNonfused[k],label=None)
                    totalRunTime += runtime
            
            elif len(filtered) and algorithm=='fused':
                foundAny=True
                
                row          = list(filter(lambda x: x[5]=='FusedTimeStep',filtered))[0]
                iterations   = float (row[10])
                runtime      = float (row[11]) / numberOfQuadraturePoints / iterations
                if i==0:
                    ax.bar(centers[i]+0.2,runtime,width=0.4,color='1.0',align='center',log=False,label='Fused')
                else:
                    ax.bar(centers[i]+0.2,runtime,width=0.4,color='1.0',align='center',log=False,label=None)
                totalRunTime = runtime
                
            yMin    = min(yMin,totalRunTime)
            yMax    = max(yMax,totalRunTime) 
    
    if foundAny:
        plt.ylabel(r'normalised time [t]=s', fontsize=8)
        plt.xlabel(r'order', fontsize=8)
        plt.grid(True, which='both')
    
        ax.yaxis.grid(False, which='minor')
        ax.xaxis.grid(False)
        ax.set_xlim([centers[0]-0.4,centers[-1]+0.4])
        ax.set_ylim([0,yMax*1.2])

        plt.tick_params(axis='both', which='major', labelsize=7)
        ax.set_xticks(centers)
        ax.set_xticklabels(list(orders)) 
        
        ax.legend(fontsize=6,loc='upper left',framealpha=0.5)   
 
        # Write files
        fig = plt.gcf()
        fig.set_size_inches(2.40,2.20) # width: 0.470 * SIAM SISC \textwidth (=5.125in)
    
        # plot
        filename = table.replace('.csv','-%s-%s-n%s-t%s-c%s' % (identifier, sharedMem, nodes, tasks, cores))
        plt.savefig('%s.pdf' % filename, bbox_inches='tight')
        plt.savefig('%s.png' % filename, bbox_inches='tight')
        print("created plot: %s.pdf" % filename)
        print("created plot: %s.png" % filename)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Compares the nonfused and fused scheme time per timestep based on the given table containing adapter 
user and cpu times (see extractadaptertimes.py)
\n\n
Sample usage:\n
python plotimestepcomparision.py -table mytable.csv 
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help='A *.csv file created with the extractadaptertimes.py script.')
parser.add_argument('-identifier',required=True,help='An identifier for the experiment. Typically the used mesh.')
parser.add_argument('-nodes',required=False,type=int,default='1',help='Number of nodes (default: 1).')
parser.add_argument('-tasks',required=False,type=int,default='1',help='Number of tasks (default: 1).')
parser.add_argument('-cores',required=False,type=int,default='1',help='Number of cores (default: 1).')
parser.add_argument('-sharedMem',required=False,type=str,default='None',help='Shared memory mode (default: None).')
parser.add_argument('-dim',required=False,type=int,default='0',help='Space dimensions. Used for normalisation (default: 0).')
parser.add_argument('-cells',required=False,type=int,default='1',help='Cells per coordinate axis. Used for normalisation (default: 1).')

args           = parser.parse_args();

table          = args.table
identifier     = args.identifier
nodes          = str(args.nodes)
tasks          = str(args.tasks)
cores          = str(args.cores)
sharedMem      = args.sharedMem
dimension      = float(args.dim)
cells          = float(args.cells)

plotTimeStepComparison(table,identifier,nodes,tasks,cores,sharedMem,dimension,cells)
