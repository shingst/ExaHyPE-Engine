#!/usr/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter
import sys

import operator
import csv

# FROM: https://stackoverflow.com/questions/17687213/how-to-obtain-the-same-font-style-size-etc-in-matplotlib-output-as-in-latex
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

def column(matrix, i):
    return [row[i] for row in matrix]

def readTable(table):
    datafile = open(table, 'r')
    next(datafile) # skip header
    data = list(csv.reader(datafile,delimiter='&'))
    datafile.close()
    return data

def query_table(data,mesh,mode,order,algorithm,cores): 
    '''
    Read a table and filter out certain
    columns 
    Args:
    table
       name of the csv file
    mesh
       mesh identifier; something like 'regular-0', 'regular-1'
    order
       approximaation order
    algorithm
       'fused' or 'nonfused'
    cores
       list of cores
    mode
        multicore mode
    '''
    return list(filter(lambda x : x[0]==mesh and x[8] == mode and int(x[1]) == order and x[4]==algorithm and int(x[7]) in cores, data))


    '''
    .. module:: usertimeplot
      :platform: Unix, Windows, Mac
      :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
       
    .. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

    :synopsis: Creates a speedup plot based on user and CPU times stored in CSV tables.
    '''
def plot_time_and_updates_per_timestep(table,mesh,mode,algorithms,cores,orders,y,ymin,ymax):
    '''
    '''
    # Plot

    for plot_option in ["runtimes", "updates"]:
        if plot_option == "updates":
            for ip in range(0,len(orders)):
                for algorithm in algorithms:
                    for i in range(0,len(y[algorithm][ip])):
                       y[algorithm][ip][i]  =  1.0 / y[algorithm][ip][i] 

        # for mode in modes:
        plt.clf() 
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        
        colors           = ['k','0.33']       # per algorithm
        markerFaceColors = ['k','None']       # per algorithm
        markers          = ['o','^','s','D']  # per order
        
        for ialg, algorithm in enumerate(algorithms):
            markerSpacing=3
            if ialg==1:
                markerSpacing=4
            
            for ip, p in enumerate(orders):
                identifier=r'$p=%d\,\rm{(%s)}$' % (p,algorithm)
                if ialg*len(orders)+ip > len(orders):
                    identifier=None
                if len(cores) == len(y[algorithm][ip]):
                    plt.plot(cores,y[algorithm][ip],lw=1.1,markersize=5,color=colors[ialg],markeredgecolor=colors[ialg],markerfacecolor=markerFaceColors[ialg],marker=markers[ip],label=identifier,markevery=markerSpacing)
                else:
                    print(algorithm,"with order",p, "is missing values")

        if plot_option == "runtimes":
            plt.ylabel(r'time per cell per real. timestep [s]', fontsize=8)
        else:
            plt.ylabel(r'cell updates per real. timestep [s]', fontsize=8)

        plt.xlabel(r'cores', fontsize=8)
        plt.grid(True, which='both')

        xticks = cores
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.axes().xaxis.set_major_locator(FixedLocator(cores[0::2]))
        plt.axes().xaxis.set_minor_locator(FixedLocator(cores[1::2]))
        
        ax.yaxis.grid(False, which='minor')
        fig = plt.gcf()
        
        # Write files
        if plot_option == "runtimes":
            ax.set_ylim(ymin*0.9,ymax*1.1)
            plt.legend(loc='upper right',fontsize=7,framealpha=0.5)
        else:
            ax.set_ylim(1.0/ymax*0.9,1.0/ymin*1.1)
            plt.legend(loc='lower right',fontsize=7,framealpha=0.5)

        fig.set_size_inches(2.40,2.20) # width: 0.470 * SIAM SISC \textwidth (=5.125in)
        
        # linear
        filename = table.replace('.runtimes.csv','-{}-{}-{}-{}'.format(mode,plot_option,mesh,'linear'))
        ax.set_yscale('linear')
        plt.savefig('%s.pdf' % filename, bbox_inches='tight')
        plt.savefig('%s.png' % filename, bbox_inches='tight')
        print('Created plot: %s.pdf' % filename)
        print('Created plot: %s.png' % filename)
        # log2
        filename = table.replace('.runtimes.csv','-{}-{}-{}-{}'.format(mode,plot_option,mesh,'log2'))
        ax.set_yscale('log',basey=2)
        plt.savefig('%s.pdf' % filename, bbox_inches='tight')
        plt.savefig('%s.png' % filename, bbox_inches='tight')
        print('Created plot: %s.pdf' % filename)
        print('Created plot: %s.png' % filename)
        # log10
        filename = table.replace('.runtimes.csv','-{}-{}-{}-{}'.format(mode,plot_option,mesh,'log10'))
        ax.set_yscale('log',basey=10)
        plt.savefig('%s.pdf' % filename, bbox_inches='tight')
        plt.savefig('%s.png' % filename, bbox_inches='tight')
        print('Created plot: %s.pdf' % filename)
        print('Created plot: %s.png' % filename)


def plot_runtimes(table,mesh,cells,dim,ignoreReruns,predictorInBackground):
    '''
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      table (str):
         A csv file
    '''
    fused    = 'fused'
    nonfused = 'nonfused'
    if predictorInBackground:
        fused    += '+bon' 
        nonfused += '+bon'
    else:     
        fused    += '+boff' 
        nonfused += '+boff'
    algorithms = [fused, nonfused]
    
    data    = readTable(table)
    orders  = sorted(list(set(map(int,column(data,1)))),key=int)
    cores   = sorted(list(set(map(int,column(data,7)))),key=int)
    modes   = sorted(list(set(column(data,8))),key=str)
    nOrders = len(orders)
    nCores  = len(cores)
    
    print(data)
    print(algorithms)
    
    normalisation = cells ** dim
    for mode in modes:
        # Read data from table and process
        ymin = sys.maxsize
        ymax = 0
        
        y = {}
        y[fused]   ={}
        y[nonfused]={}
        
        foundOrders = orders
        for ip, p in enumerate(orders):
            dataFused    = query_table(data,mesh,mode,p,fused,cores)
            dataNonfused = query_table(data,mesh,mode,p,nonfused,cores)
            
            if len(dataFused) and len(dataNonfused):
                y[nonfused][ip]  = list( map( float, column(dataNonfused, 10) ) )
                if ignoreReruns is True:
                    y[fused][ip] = list( map( float, column(dataFused, 11) ) )
                else:
                    y[fused][ip] = list( map( float, column(dataFused, 10) ) )
                # Normalise fused values w.r.t. number of cells
                itFused    = float(dataFused[0][9])
                itNonfused = float(dataNonfused[0][9])
                for i in range(0,len(y[fused][ip])):
                    y[fused][ip][i] /= normalisation
                    y[fused][ip][i] *= itFused / itNonfused # Normalise fused values wrt iterations of nonfused
                for i in range(0,len(y[nonfused][ip])):
                    y[nonfused][ip][i] /=normalisation
                
                for algorithm in algorithms:
                    ymax = max( ymax, max(y[algorithm][ip]) )
                    ymin = min( ymin, min(y[algorithm][ip]) )
            else:
                foundOrders.remove(p)
        
        # Plot
        plot_time_and_updates_per_timestep(table,mesh,mode,algorithms,cores,foundOrders,y,ymin,ymax) 

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
parser.add_argument('-mesh',required=True,help='Mesh identifier.')
parser.add_argument('-cells',required=True,help='Number of cells per axis')
parser.add_argument('-dim',required=True,help='Space dimensions.')
parser.add_argument('--ignore-reruns',dest='ignoreReruns',required=False,action='store_true',help="Do not consider predictor reruns in the runtime calculations.")
parser.add_argument('--predictor-in-background',dest='predictorInBackground',required=False,action='store_true',help="Predictor is run as background thread.")
parser.set_defaults(ignoreReruns=False)
parser.set_defaults(predictorInBackground=False)

args = parser.parse_args();

table                 = args.table
mesh                  = args.mesh
cells                 = int(args.cells)
dim                   = int(args.dim)
ignoreReruns          = args.ignoreReruns
predictorInBackground = args.predictorInBackground

plot_runtimes(table,mesh,cells,dim,ignoreReruns,predictorInBackground)
