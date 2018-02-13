#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

'''
.. module:: extractadaptertimes
  :platform: Unix, Windows, Mac
  :synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
'''
def parse_adapter_times(file_path,per_iteration=False):
    '''
    Reads a single Peano output file and parses the user time spent within each adapter.
    
    Args:
       file_path (str):
          Path to the Peano output file.
       per_iteration (bool):
          Specifies if the times per iteration should be read in instead of the total times.

    Returns:
       A dict holding for each of the found adapters a nested dict that holds the following key-value pairs:
          * 'n'       : (int)    Number of times this adapter was used.
          * 'cputime' : (float) CPU time spent within the adapter.
          * 'usertime': (float) User time spent within the adapter.
    '''
    result = { }
    try:
        file_handle=open(file_path)
        
        cputime_index  = 3
        usertime_index = 5
        if per_iteration:
            cputime_index  = 4
            usertime_index = 6
        
        for line in file_handle:
            anchor = '|'
            header = '||'
            
            if anchor in line and header not in line:
                segments = line.split('|')
                adapter = segments[1].strip();
                result[adapter]             = {}
                result[adapter]['n']        = int(segments[2].strip())
                result[adapter]['cputime']  = float(segments[cputime_index ].strip())
                result[adapter]['usertime'] = float(segments[usertime_index].strip())
    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise
    return result

def extract_singlecore_table(root_dir,prefix):
    '''
    Extracts adapter times from Peano output files with specific file naming pattern.
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
    header = ["Order","Adapter","Optimisation","Run","Iterations","User Time (Total)","CPU Time (Total)"]

    # collect filenames
    with open(root_dir+"/"+prefix+'.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        # write header
        csvwriter.writerow(header)
        
        # write content
        print("processed files:")
        for filename in os.listdir(root_dir):
            if filename.endswith(".out") and filename.startswith(prefix):
                # file structure: (prefix)-(arch)-(optimisation)-p(order)-r(run).out
                match = re.search('^'+prefix+'-(.+)-p([0-9]+)-r([0-9]+)\.out',filename)
                print(root_dir+"/"+filename)
                
                optimisation = match.group(1) # opt/gen
                order        = match.group(2) 
                run          = match.group(3) 
                    
                times = parse_adapter_times(root_dir+"/"+filename) 
                
                for adapter in times:
                    iterations = times[adapter]['n']
                    usertime   = times[adapter]['usertime']
                    cputime    = times[adapter]['cputime']

                    csvwriter.writerow([order,adapter,optimisation,run,iterations,usertime,cputime]) 

def extract_table(root_dir,prefix):
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
    header = ["Mesh","Order","CC","Kernels","Algorithm","Adapter","Nodes","Tasks (per Node)","Cores (per Task)","Mode","Iterations","User Time (Total)","CPU Time (Total)"]

    # collect filenames
    with open(root_dir+"/"+prefix+'.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        # write header
        csvwriter.writerow(header)
        
        # write content
        print("processed files:")
        for filename in os.listdir(root_dir):
            if filename.endswith(".out") and filename.startswith(prefix):
                # sample: Elastic3D-no-output-gen-fused-regular-0-p3-TBB-Intel-n1-t1-c24.out
                match = re.search('^'+prefix+'-([a-z]+)-(([A-Za-z]|\+)+)-(.+)-p([0-9]+)-(.+)-([A-Za-z]+)-n([0-9]+)-t([0-9]+)-c([0-9]+).out$',filename)
                print(root_dir+"/"+filename)
                kernels   = match.group(1) # opt/gen
                algorithm = match.group(2) # fused/nonfused
                mesh      = match.group(4)
                order     = match.group(5)
                mode      = match.group(6)
                cc        = match.group(7)
                nodes     = match.group(8)
                tasks     = match.group(9)
                cores     = match.group(10)
                    
                times = parse_adapter_times(root_dir+"/"+filename) 
                
                for adapter in times:
                    iterations = times[adapter]['n']
                    usertime   = times[adapter]['usertime']
                    cputime    = times[adapter]['cputime']

                    row=[mesh,order,cc,kernels,algorithm,adapter,nodes,tasks,cores,mode,iterations,usertime,cputime]
                    csvwriter.writerow(row)

def sort_table(filename,myKey):
    '''
    Sorts the rows of the file according to nodes,tasks,cores,adapter name.
    See: https://stackoverflow.com/a/17109098
    '''
    datafile    = open(filename, 'r')
    header      = next(datafile).strip()
    reader      = csv.reader(datafile,delimiter='&')
    # [order,cc,kernels,algorithm,adapter,nodes,tasks,cores,mode,iterations,usertime,cputime]
    sorted_data = sorted(reader, key=myKey)
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
<prefix>.csv

\n\n
Sample usage:\n
python3 extractadaptertimes.py -path \'results/' -prefix \'Euler-ADERDG-regular-0-fused-p3\'
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the Peano output files.")
parser.add_argument('--single-core',dest='singlecore',required=False,action='store_true',help="")
parser.set_defaults(singlecore=False)

args     = parser.parse_args();

root_dir   = args.path
prefix     = args.prefix
singlecore = args.singlecore 

if singlecore is False:
  extract_table(root_dir,prefix)
  key = lambda x: (x[0],int(x[1]),x[2],x[3],x[4],x[5],int(x[6]),int(x[7]),int(x[8])) 
  sort_table(root_dir+"/"+prefix+".csv",key)
else:
  extract_singlecore_table(root_dir,prefix)
  key = lambda x: ( int(x[0]), x[1], x[2], x[3], int(x[4]) ) 
  sort_table(root_dir+"/"+prefix+".csv",key)

print("created table:")
print(root_dir+"/"+prefix+".csv")
