#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

'''
.. module:: extractlikwidmetrics
  :platform: Unix, Windows, Mac
  :synopsis: Contains routines to extract 
             performance metrics from Peano output files with specific file naming pattern.
             The call of the Peano application must be wrapped within likwid-perfctr via: 
             likwid-perfctr -C 0 -g <METRIC> <application>
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
'''

def extract_likwid_metrics(root_dir,prefix):
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
    
    We are currently parsing  "MFLOP/s STAT" (SUM), "Memory bandwidth [MBytes/s]" (SUM),
    "Branch misprediction rate STAT" (Avg),  "L2 miss rate STAT" (Avg).
    
    The output of likwid is different for singlecore and multicore
    runs:
    
    # Singlecore Metrics

    ## FLOPS_DP

    +----------------------+-----------+
    |        Metric        |   Core 1  |
    +----------------------+-----------+
    |  Runtime (RDTSC) [s] |   6.4970  |
    | Runtime unhalted [s] |   7.7641  |
    |      Clock [MHz]     | 2867.1614 |
    |          CPI         |   1.4365  |
    |        MFLOP/s       |  578.7295 |
    |      AVX MFLOP/s     |     0     |
    |    Packed MUOPS/s    |  252.4243 |
    |    Scalar MUOPS/s    |  73.8808  |
    +----------------------+-----------+

    ##  MEM

    +-----------------------------------+------------+
    |               Metric              |   Core 1   |
    +-----------------------------------+------------+
    |        Runtime (RDTSC) [s]        |   6.5468   |
    |        Runtime unhalted [s]       |   7.7608   |
    |            Clock [MHz]            |  2849.4441 |
    |                CPI                |   1.4359   |
    |  Memory read bandwidth [MBytes/s] |  9203.1065 |
    |  Memory read data volume [GBytes] |   60.2505  |
    | Memory write bandwidth [MBytes/s] |  4409.1727 |
    | Memory write data volume [GBytes] |   28.8658  |
    |    Memory bandwidth [MBytes/s]    | 13612.2793 |
    |    Memory data volume [GBytes]    |   89.1163  |
    +-----------------------------------+------------+

    ## BRANCH

    +----------------------------+--------------+
    |           Metric           |    Core 1    |
    +----------------------------+--------------+
    |     Runtime (RDTSC) [s]    |    6.5219    |
    |    Runtime unhalted [s]    |    7.7536    |
    |         Clock [MHz]        |   2858.6786  |
    |             CPI            |    1.4345    |
    |         Branch rate        |    0.1190    |
    |  Branch misprediction rate | 2.575205e-07 |
    | Branch misprediction ratio | 2.163851e-06 |
    |   Instructions per branch  |    8.4026    |
    +----------------------------+--------------+

    ## L2CACHE

    +----------------------+-----------+
    |        Metric        |   Core 1  |
    +----------------------+-----------+
    |  Runtime (RDTSC) [s] |   6.5188  |
    | Runtime unhalted [s] |   7.7477  |
    |      Clock [MHz]     | 2855.7219 |
    |          CPI         |   1.4335  |
    |    L2 request rate   |   0.3787  |
    |     L2 miss rate     |   0.1056  |
    |     L2 miss ratio    |   0.2788  |
    +----------------------+-----------+


    # Multicore Metrics


    ## FLOPS_DP:

    +---------------------------+------------+-----------+-----------+-----------+
    |           Metric          |     Sum    |    Min    |    Max    |    Avg    |
    +---------------------------+------------+-----------+-----------+-----------+
    |  Runtime (RDTSC) [s] STAT |   27.7776  |   1.1574  |   1.1574  |   1.1574  |
    | Runtime unhalted [s] STAT |   17.5755  |   0.7025  |   1.1664  |   0.7323  |
    |      Clock [MHz] STAT     | 53287.5604 | 2194.3481 | 2339.6363 | 2220.3150 |
    |          CPI STAT         |   80.7068  |   0.8296  |   3.5650  |   3.3628  |
    |        MFLOP/s STAT       |  3248.7237 |  118.0829 |  532.8151 |  135.3635 |
    |      AVX MFLOP/s STAT     |      0     |     0     |     0     |     0     |
    |    Packed MUOPS/s STAT    |  1416.9944 |  59.0414  |  59.0415  |  59.0414  |
    |    Scalar MUOPS/s STAT    |  414.7347  |     0     |  414.7321 |  17.2806  |
    +---------------------------+------------+-----------+-----------+-----------+

    ## MEM:

    +----------------------------------------+------------+-----------+------------+-----------+
    |                 Metric                 |     Sum    |    Min    |     Max    |    Avg    |
    +----------------------------------------+------------+-----------+------------+-----------+
    |        Runtime (RDTSC) [s] STAT        |   28.1040  |   1.1710  |   1.1710   |   1.1710  |
    |        Runtime unhalted [s] STAT       |   17.6412  |   0.6944  |   1.1724   |   0.7351  |
    |            Clock [MHz] STAT            | 53428.9064 | 2153.7694 |  2358.9159 | 2226.2044 |
    |                CPI STAT                |   81.0719  |   0.8342  |   3.5820   |   3.3780  |
    |  Memory read bandwidth [MBytes/s] STAT | 52486.4902 |     0     | 26297.3151 | 2186.9371 |
    |  Memory read data volume [GBytes] STAT |   61.4616  |     0     |   30.7941  |   2.5609  |
    | Memory write bandwidth [MBytes/s] STAT | 26486.8194 |     0     | 13702.7744 | 1103.6175 |
    | Memory write data volume [GBytes] STAT |   31.0160  |     0     |   16.0459  |   1.2923  |
    |    Memory bandwidth [MBytes/s] STAT    | 78973.3098 |     0     | 40000.0896 | 3290.5546 |
    |    Memory data volume [GBytes] STAT    |   92.4777  |     0     |   46.8401  |   3.8532  |
    +----------------------------------------+------------+-----------+------------+-----------+

    ## BRANCH:

    +---------------------------------+------------+--------------+--------------+--------------+
    |              Metric             |     Sum    |      Min     |      Max     |      Avg     |
    +---------------------------------+------------+--------------+--------------+--------------+
    |     Runtime (RDTSC) [s] STAT    |   27.2616  |    1.1359    |    1.1359    |    1.1359    |
    |    Runtime unhalted [s] STAT    |   17.7047  |    0.7118    |    1.1711    |    0.7377    |
    |         Clock [MHz] STAT        | 53381.5227 |   2197.7583  |   2378.3713  |   2224.2301  |
    |             CPI STAT            |   81.2767  |    0.8336    |    3.5603    |    3.3865    |
    |         Branch rate STAT        |   3.0690   |    0.0450    |    0.1340    |    0.1279    |
    |  Branch misprediction rate STAT |   0.0001   | 9.050274e-07 | 2.820245e-05 | 5.374324e-06 |
    | Branch misprediction ratio STAT |   0.0010   | 6.977308e-06 |    0.0002    | 4.115248e-05 |
    |   Instructions per branch STAT  |  197.2038  |    7.4615    |    22.2424   |    8.2168    |
    +---------------------------------+------------+--------------+--------------+--------------+

    ## L2CACHE

    +---------------------------+------------+-----------+-----------+-----------+
    |           Metric          |     Sum    |    Min    |    Max    |    Avg    |
    +---------------------------+------------+-----------+-----------+-----------+
    |  Runtime (RDTSC) [s] STAT |   27.4632  |   1.1443  |   1.1443  |   1.1443  |
    | Runtime unhalted [s] STAT |   17.7995  |   0.7182  |   1.1748  |   0.7416  |
    |      Clock [MHz] STAT     | 53545.8000 | 2203.7171 | 2357.4567 | 2231.0750 |
    |          CPI STAT         |   81.9969  |   0.8364  |   3.5848  |   3.4165  |
    |    L2 request rate STAT   |   12.3080  |   0.1178  |   0.5385  |   0.5128  |
    |     L2 miss rate STAT     |   2.9400   |   0.0293  |   0.1286  |   0.1225  |
    |     L2 miss ratio STAT    |   5.7406   |   0.2358  |   0.2489  |   0.2392  |
    +---------------------------+------------+-----------+-----------+-----------+
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
    
    metrics    = [
                  ["  MFLOP/s",                   "Sum"],   # Two whitespaces are required to not find the AVX MFLOP/s by accident
                  ["AVX MFLOP/s",                 "Sum"],
                  ["Memory bandwidth [MBytes/s]", "Sum"],
                  ["Memory data volume [GBytes]", "Sum"],
                  ["L3 bandwidth [MBytes/s]",     "Sum"], 
                  ["L3 data volume [GBytes]",     "Sum"],
                  ["L3 request rate",             "Avg"],
                  ["L3 miss rate",                "Avg"],
                  ["L2 request rate",             "Avg"],
                  ["L2 miss rate",                "Avg"],
                  ["Branch misprediction rate",   "Avg"]
                 ]
    
    # collect filenames
    with open(root_dir+"/"+prefix+'.likwid.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        # write header
        header = ["Mesh","Order","CC","Kernels","Algorithm","Nodes","Tasks (per Node)","Cores (per Task)","Shared Memory"]
        for metric in metrics:
            header.append(metric[0]+"("+metric[1]+")")
        csvwriter.writerow(header)

        # write content
        print("processed files:")
        for filename in os.listdir(root_dir):
            if filename.endswith(".out.likwid") and filename.startswith(prefix):
                print(root_dir+"/"+filename)
                # sample: Elastic3D-no-output-gen-fused-regular-0-p3-TBB-Intel-n1-t1-c24.out
                match = re.search('^'+prefix+'-([a-z]+)-([a-z]+)-(.*)-p([0-9]+)-([A-Za-z]+)-([A-Za-z]+)-n([0-9]+)-t([0-9]+)-c([0-9]+)',filename)
                kernels   = match.group(1) # opt/gen
                algorithm = match.group(2) # fused/nonfused
                mesh      = match.group(3)
                order     = match.group(4)
                mode      = match.group(5)
                cc        = match.group(6)
                nodes     = match.group(7)
                tasks     = match.group(8)
                cores     = match.group(9)
                    
                measurements = parse_likwid_metrics(root_dir+'/'+filename,metrics,int(cores)==1) 
                
                row = [mesh,order,cc,kernels,algorithm,nodes,tasks,cores,mode]
                   
                for metric in metrics:
                    row.append ( str(measurements[metric[0]][metric[1]]) )
                csvwriter.writerow(row)

def parse_likwid_metrics(file_path,metrics,singlecore=False):
    """
    Reads a single Peano output file and parses likwid performance metrics.
    
    Args:
       file_path (str):
          Path to the Peano output file.
       metrics (str[][]):
          A list of metrics the we want to read out.
       singlecore (bool):
          Specifies if the run was a singlecore run.

    Returns:
       A dict holding for each of the found metrics a nested dict that holds the following key-value pairs:
          * 'Sum' 
          * 'Avg' 
          * 'Min' 
          * 'Max' 
    """
    columns    = [ "Sum","Min","Max","Avg" ]
    
    result  = { }
    for metric in metrics:
        result[metric[0]] =  { }
        result[metric[0]][metric[1]] = -1.0

    try:
        file_handle=open(file_path)
        
        for line in file_handle:
            for metric in metrics: 
                if singlecore:
                    if metric[0] in line:
                        segments = line.split('|')
                        
                        #    |     Runtime (RDTSC) [s]    |    6.5219    |
                        value  = float(segments[2].strip());
                        values = {}                         
                        values["Sum"] = value
                        values["Min"] = value
                        values["Max"] = value
                        values["Avg"] = value
                        result[metric[0]][metric[1]]=values[metric[1]]                        
                else:
                    if metric[0]+" STAT" in line:
                        segments = line.split('|')
                        #   |  Runtime (RDTSC) [s] STAT |   27.4632  |   1.1443  |   1.1443  |   1.1443  |
                        values = {}                                                 
                        values["Sum"] = float(segments[2].strip());
                        values["Min"] = float(segments[3].strip());
                        values["Max"] = float(segments[4].strip());
                        values["Avg"] = float(segments[5].strip());
                        result[metric[0]][metric[1]]=values[metric[1]]
    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise
    return result


def sort_table(filename):
    '''
    Sorts the rows of the file according to nodes,tasks,cores,
    See: https://stackoverflow.com/a/17109098
    '''
    datafile    = open(filename, 'r')
    header      = next(datafile).strip()
    reader      = csv.reader(datafile,delimiter='&')
    # row = [order,cc,kernels,algorithm,nodes,tasks,cores,mode]
    sorted_data = sorted(reader, key=lambda x: (x[0],int(x[1]),x[2],x[3],x[4],int(x[5]),int(x[6]),int(x[7])))
    datafile.close() 
 
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header.split('&'))
        writer.writerows(sorted_data)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract performance metrics from Peano output files with specific file naming pattern
and write them to a csv file with name
<prefix.csv

\n\n
Sample usage:\n
python extractlikwidmetrics.py -path \'examples/151217_phi1_node/' -prefix "Euler"
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the Peano output files.")

args     = parser.parse_args();

root_dir = args.path
prefix   = args.prefix

extract_likwid_metrics(root_dir,prefix)
sort_table(root_dir+"/"+prefix+".likwid.csv")
print("created table:")
print(root_dir+"/"+prefix+".likwid.csv")
