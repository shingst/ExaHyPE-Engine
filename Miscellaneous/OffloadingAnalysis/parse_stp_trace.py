#!/usr/bin/env python3

import sys
import os

#basename = "exahype_solvers_ADERDGSolver_MigratablePredictionJob__run_iterations_rank_"

def readFile(filename):
  accIterations = 0
  accTime = 0
  
  if os.path.isfile(filename):

    f = open(filename, "r")
    for line in f:
      splitted = (line[:-1].split(":"))
      run = int(splitted[0])
      iterations = int(splitted[1])
      accIterations += iterations+1
      accTime += run
    f.close()

  return (accIterations, accTime)

def parseDir(dirname, ranks, threads, timestep, basename):
  totaliterations = []
  totaltime = []

  x = range(1,ranks)

  for i in range(1,ranks): 
    accIterationsRank = 0
    accTimeRank = 0

    for j in range(0, threads):
      filename =  dirname+"/"+basename+str(i)+"_"+str(j)+"_step_"+str(timestep)+".txt"
      print ("processing", filename)
      (its, run) = readFile(filename)
      #print(its, run)
      accIterationsRank += its
      accTimeRank += run

    totaliterations.append(accIterationsRank)
    totaltime.append(accTimeRank)
  return (x, totaliterations, totaltime)
  

if __name__=="__main__":
  dirname = sys.argv[1]
  ranks = int(sys.argv[2])
  threads = int(sys.argv[3])


  (x, totaliterations, totaltime) = parseDir(directory, ranks, threads)

  print (x)
  print (totaliterations)
  print (totaltime)
