#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re

'''
.. module:: extracterrors
  :platform: Unix, Windows, Mac
  :synopsis: Extracts errors from Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts errors from Peano output files with specific file naming pattern.
'''


def interpolate_errors(file_path,t_eval):
    """
    Reads a single Peano output file, reads the written out errors and interpolates
    the errors linearly in order to determine  t_eval.
    
    Args:
       file_path (str):
          Path to the Peano output file.
       t_eval (float):
          Evaluation time
    """
    errors = ["absErrorL1","absErrorL2","absErrorLInf","relErrorL1","relErrorL2","relErrorLInf"]

    variables=5

    # init result set   
    try:
        file_handle=open(file_path)
       
        hMax    = 0. 
        t1      = 0. 
        t2      = 0.
        
        result1 = {}        
        for variable in range(0,variables):
            result1[variable] = {}
            for error in errors:
                result1[variable][error] = 0.0
        
        result2 = {}        
        for variable in range(0,variables):
            result2[variable] = {}
            for error in errors:
                result2[variable][error] = 0.0

        t=0
        for line in file_handle:
            if "coarsest mesh size was chosen as " in line:
                match = re.search('chosen as ([0-9]+\.[0-9]+)',line)
                hMax = float(match.group(1).strip())
            
            if "t_eval" in line:
                t = float(line.split(':')[-1].strip())
                if t2>0:
                    dt = t2-t1
                
                    print ("t_eval=%1.6g"%t_eval) 
                    print("t_before=%1.6g"%t1)
                    print("t_after=%1.6g"%t2)
                    if dt>0 and t1>0:
                        for error in errors:
                            for variable in range(0,variables):
                                interpoland = result1[variable][error] + ( result2[variable][error]-result1[variable][error] ) / dt * (t_eval-t1)
                                print (error+"(%d)= %1.2e" % (variable, interpoland))
                    else: 
                        for error in errors:
                            for variable in range(0,variables):
                                interpoland = result2[variable][error]
                                print (error+"(%d)= %1.2e" % (variable, interpoland))
                    break # line in file_handle

            # interpolate            
            if t <= t_eval:
                t1 = t
                # parser errors
                for error in errors:
                    if error in line:
                        segments = line.split(':')[1].split(',')
                    
                        for variable in range(0,variables):
                            result1[variable][error] = float(segments[variable])
                 

            if t >= t_eval: 
                t2 = t
                # parser errors
                for error in errors:
                    if error in line:
                        segments = line.split(':')[1].split(',')
                    
                        for variable in range(0,variables):
                            result2[variable][error] = float(segments[variable])
             

    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Reads a single Peano output file, reads the written out errors and interpolates
the errors linearly in order to determine  t_eval.
\n\n
Sample usage:\n
python3 extracterrors.py -path \'results/' -prefix \'Euler-ADERDG-regular-0-fused-p3\'
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="path to a file")
parser.add_argument('-t',required=True,help="Evaluation time.")

args = parser.parse_args();

path   = args.path
t_eval = args.t

interpolate_errors(path,float(t_eval))
