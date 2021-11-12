#!/bin/python3
import os
import sys
import re

pattern_success = re.compile(".*Peano terminates successfully.*")
pattern_soft_error_detected = re.compile(".*soft error detected.*")

def main():
 dir = sys.argv[1]
 runs = int(sys.argv[2])

 success = 0
 detected = 0

 success_run = 0
 detected_run = 0

 for r in range(1,runs+1):
   for filename in os.listdir(dir):
     if filename.endswith("-r"+str(r)+".out"):
       for i, line in enumerate(open(os.path.join(dir,filename))):
         if re.match(pattern_success, line):
           success_run += 1
     elif filename.endswith("-r"+str(r)+".job_err"):
       for i, line in enumerate(open(os.path.join(dir,filename))):
         if re.match(pattern_soft_error_detected, line):
           detected_run += 1
           break

   print("run ", r, success_run, detected_run)

   detected += detected_run
   success += success_run
   
   success_run = 0
   detected_run = 0


 print("successful runs:", success/2)
 print("detected errors:", detected)

if __name__ == "__main__":

  main()
