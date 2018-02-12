#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import sys
import re
import os
import csv

import operator

def column(matrix, i):
    return [row[i] for row in matrix]

def readTable(filename):
    datafile = open(filename, 'r')
    next(datafile) # skip header
    data = list(csv.reader(datafile,delimiter='&'))
    datafile.close()
    return data

def writeTable(data,header,filename):
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        
        for row in data:
            writer.writerow(row)
    
def extractSingleCoreRuntimes(table):
    '''
    assume: inputRow = [order,adapter,optimisation,run,iteratins,user time, runtime]
    example:
       0           1             2            3     4        5         6
    ['6', 'MeshRefinement', 'noarch-O3-vec', '1', '12', '15.2857', '15.25'] 
    
    Args:
      table (str):
         path of a CSV table created with extractadaptertimes.py
    '''
    nonfusedAdapters = [ 
                         "BroadcastGlobalDataAndMergeNeighbourMessages",
                         "SolutionUpdate",
                         "Prediction" 
                       ]
    data = readTable(table)
    
    orders        = set(column(data,0))
    optimisations = set(column(data,2))
    
    # processing data
    result = []
    for order in orders:
        for optimisation in optimisations:
            averageRuntimePerIteration = 0.0
            minIterations              = sys.maxsize
            for adapter in nonfusedAdapters:
                filtered = list(filter(lambda x: x[0]==order and x[1]==adapter and x[2]==optimisation, data))
                runs     = len(filtered)
                
                iterations = 0
                if runs > 0:
                    iterations = int ( filtered[0][4] )
                    
                    averageAdapterTime = 0.0
                    for r in range(0,runs):
                        averageAdapterTime += float ( filtered[r][5] )
                    averageAdapterTime /= runs
                    
                    averageRuntimePerIteration += averageAdapterTime / iterations
                minIterations = min(minIterations,iterations)
                
            row = [order,optimisation,minIterations,averageRuntimePerIteration]
            result.append(row)
    
    header = ["Order","Optimisation","Iterations","Total Runtime (Per Iteration)"]
    filename = table.replace('.csv','.runtimes.csv')
    result = sorted(result, key=lambda x: (int(x[0]),x[1]))
    writeTable(result,header,filename)
'''
.. module:: extractunigridruntimes
  :platform: Unix, Windows, Mac
  :synopsis: Extract accumulated adapter runtimes from CSV files created
             by the extractadaptertimes.py script.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extract accumulated adapter runtimes from CSV files created
           by the extractadaptertimes.py script.
'''

def extractRuntimes(table):
    '''
    Args:
      table (str):
         path of a CSV table created with extractadaptertimes.py
      single-core (bool)
         If true, use old adapter names: ADERDGTimeStep,... Not the new ones: FusedTimeStep,...
    '''
    nonfusedAdapters = [ 
                          'BroadcastGlobalDataAndMergeNeighbourMessages',
                          'SolutionUpdate',
                          'Prediction' # substract one; exclude initialisation
                       ]
    fusedAdapters    = [ 
                        'FusedTimeStep',
                        'PredictionRerun'
                       ]
                       
    
    data = readTable(table)
    
    #                     0       1    2   3         4       5       6     7     8     9    10            11              12 
    # assume: row = [identifier,order,cc,kernels,algorithm,adapter,nodes,tasks,cores,mode,iterations,total user time, total cpu time]
    data = sorted(data, key=lambda x: (x[0],x[1],x[2],x[3],x[4],int(x[6]),int(x[7]),int(x[8]),x[5]))
    
    identifiers = sorted(list(set(column(data,0))),key=str)
    orders      = sorted(list(set(column(data,1))),key=int)
    compilers   = sorted(list(set(column(data,2))),key=str)
    kernels     = sorted(list(set(column(data,3))),key=str)
    algorithms  = sorted(list(set(column(data,4))),key=str)
    nodes       = sorted(list(set(column(data,6))),key=int)
    tasks       = sorted(list(set(column(data,7))),key=int)
    cores       = sorted(list(set(column(data,8))),key=int)
    modes       = sorted(list(set(column(data,9))),key=str)

    
    # processing data
    result = []
    for identifier in identifiers:
        for order in orders:
            for compiler in compilers:
                for kernel in kernels:
                    for node in nodes:
                        for task in tasks:
                            for core in cores:
                                for mode in modes:
                                    for algorithm in algorithms:
                                        averageRuntimePerIteration       = 0.0
                                        averageRuntimePerFusedTimeStep   = 0.0
                                        averageRuntimePerPredictionRerun = 0.0
                                        
                                        timesteps        = 0
                                        predictionReruns = 0
                                        found = False
                                        if algorithm.startswith('fused'):
                                            for adapter in fusedAdapters:
                                                filtered = list(filter(lambda x: \
                                                    x[0]==identifier and x[1]==order and x[2]==compiler and x[3]==kernel and \
                                                    x[4]==algorithm and x[5]==adapter and x[6]==node and x[7]==task and \
                                                    x[8]==core and x[9] == mode, data))
                                                
                                                if len(filtered):
                                                    found = True
                                                    iterations = int ( filtered[0][10] )
                                                    if adapter=='FusedTimeStep':
                                                        timesteps                        = iterations
                                                        averageRuntimePerFusedTimeStep   = float ( filtered[0][11] )
                                                    elif adapter=='PredictionRerun':
                                                        predictionReruns                = iterations
                                                        averageRuntimePerPredictionRerun = float ( filtered[0][11] )

                                            if timesteps > 0:
                                                averageRuntimePerIteration        = (averageRuntimePerFusedTimeStep+averageRuntimePerPredictionRerun) / timesteps
                                                averageRuntimePerFusedTimeStep   /= timesteps
                                            if predictionReruns > 0:
                                                averageRuntimePerPredictionRerun /= predictionReruns                    
                             
                                        
                                        elif algorithm.startswith('nonfused'):
                                            averageRuntimePerIteration = 0.0
                                            timesteps = sys.maxsize
                                            for adapter in nonfusedAdapters:
                                                filtered = list(filter(lambda x: \
                                                    x[0]==identifier and x[1]==order and x[2]==compiler and x[3]==kernel and \
                                                    x[4]==algorithm and x[5]==adapter and x[6]==node and x[7]==task and \
                                                    x[8]==core and x[9] == mode, data))

                                                if len(filtered):
                                                    found = True
                                                    iterations                  = int    ( filtered[0][10] )
                                                    averageRuntimePerIteration += float ( filtered[0][11] ) / iterations
                                                    timesteps                   = min(timesteps,iterations)
                                         
                                        if found:        
                                            row = [identifier,order,compiler,kernel,algorithm,node,task,core,mode,timesteps,averageRuntimePerIteration,averageRuntimePerFusedTimeStep,predictionReruns,averageRuntimePerPredictionRerun]
                                            result.append(row)
    
    # loop
    header = ["Mesh","Order","CC","Kernels","Algorithm","Nodes","Tasks (per Node)","Cores (per Task)","Mode","Iterations","Total Runtime (Per Iteration)","FusedTimeStep (Per Iteration)","Reruns", "Predictor Rerun (Per Iteration)"]
    filename = table.replace('.csv','.runtimes.csv')
    # print(result[0])
    result = sorted(result, key=lambda x: (x[0],int(x[1]),x[2],x[3],x[4],int(x[5]),int(x[6]),int(x[7])))
    writeTable(result,header,filename)
########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract accumulated adapter runtimes from CSV tables created 
by the extractadaptertimes.py script.
Write the result to another CSV file with name
<table>.runtimes.csv

NOTE: 
This script is solely useful for computing
the runtimes of the fused and nonfused ADER-DG
algorithms on uniform meshes.
Does ignore the mesh setup times.

NOTE: Assumes no plotting was performed during
the simulations!

NOTE: Does not consider AMR or Limiting.
Parsed solve times could be completely wrong!

NOTE: Does not sanitise the input. Be careful that
there are no whitelines in the table file!

\n\n
Sample usage:\n
python extractruntimes.py -table Elastic3D-no-output.csv --single-core
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help="Directory containing the Peano output files.")
parser.add_argument('--single-core',dest='singlecore',required=False,action='store_true',help="Read a single-core table which is structured differently than the other tables.")
parser.set_defaults(singlecore=False)

args = parser.parse_args();

table      = args.table
singlecore = args.singlecore

if singlecore is False:
  extractRuntimes(table)
else:
  extractSingleCoreRuntimes(table)  

print("created table:")
print(table.replace('.csv','.runtimes.csv'))
