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
  :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on user and CPU times stored in CSV tables.
'''

def column(matrix, i):
    return [row[i] for row in matrix]

def plotTimePerTimeStep(table,dimension,cells):
    '''
    '''
    # Plot
    fig = plt.figure()
    ax  = fig.add_subplot(111)

    datafile = open(table, 'r')
    next(datafile) # skip header
    data = list(csv.reader(datafile,delimiter='&'))
    datafile.close()
    
    orders        = sorted(list(set(column(data,0))),key=str)
    optimisations = sorted(list(set(column(data,1))),key=str)
    
    print(optimisations)
    
    centers = range(0,len(orders))
    width   = 0.9 / len(optimisations)
    
    N    = len(optimisations)
    yMin = 10.0**20
    yMax = 0.0
    
    for i,order in enumerate(orders):
        DOFS = ( cells*(int(order)+1)**dimension ) # TIME PER Q
        
        for j,optimisation in enumerate(optimisations):
            row = list(filter(lambda x: x [0]==order and x[1]==optimisation, data))
            
            runtime = float(row[0][3]) / DOFS
            yMin    = min(yMin,runtime)
            yMax    = max(yMax,runtime)
            
            colour = ( len(optimisations)-1 -float(j) ) / (len(optimisations)-1)
            ax.bar(centers[i]-width*N/2+j*width, runtime,width=width,color=str(colour),align='center',log=False)
        
    plt.ylabel(r'normalised time [t]=s', fontsize=8)
    plt.xlabel(r'order', fontsize=8)
    plt.grid(True, which='both')
    
    plt.tick_params(axis='both', which='major', labelsize=7)
    
    ax.set_xticks(centers)
    ax.set_xticklabels(list(orders))
    
    ax.yaxis.grid(False, which='minor')
    ax.xaxis.grid(False)
    ax.set_xlim(-(N+1)/2*width,len(orders)-1+(N+1)/2*width)
    ax.set_ylim([yMin*0.9,yMax*1.1])
    
    # Write files
    fig = plt.gcf()
    fig.set_size_inches(2.40,2.20) # width: 0.470 * SIAM SISC \textwidth (=5.125in)
    
    # plot
    filename = table.replace('.runtimes.csv','-singlecore-runtimes')
    plt.savefig('%s.pdf' % filename, bbox_inches='tight')
    plt.savefig('%s.png' % filename, bbox_inches='tight')

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on the given tables containing adapter user and cpu times.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n
python plotmulticorespeedup.py -table 
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',required=True,help='A .runtimes.csv file created with the extractruntimes.py script.')
parser.add_argument('-dimension',required=True,help='Space dimensions')
parser.add_argument('-cells',required=True,help='Cells per coordinate axis')

args           = parser.parse_args();

table          = args.table
dimension      = float(args.dimension)
cells          = float(args.cells)

plotTimePerTimeStep(table,dimension,cells)
