#!/usr/bin/env python3
"""
.. module:: sweep_analysis
  :platform: Unix, Windows, Mac
  :synopsis: Submodule containing modules for analysing a bunch of
   ExaHyPE output files
   and creating

.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>,

:synopsis: Generate benchmark suites for ExaHyPE.
"""
import csv
import json
import sys,os
import re
import codecs
import math

import collections
import statistics

import argparse

lastKeyColumn=-1

knownParameters  = ["architecture", "dimension"]

# |      DP MFLOP/s STAT      |   9041.3374 |   68.6500 |  415.9467 |  251.1483 |
# |    AVX DP MFLOP/s STAT    |   8052.8481 |   61.4873 |  377.0757 |  223.6902 |
# |   AVX512 DP MFLOP/s STAT  |   5974.9405 |   45.5262 |  279.9348 |  165.9706 |

metrics =  [
        ["  MFLOP/s",                   "Sum"],  # Two whitespaces are required to not find the AVX MFLOP/s by accident
        ["AVX MFLOP/s",                 "Sum"],
        ["  DP MFLOP/s",                "Sum"],  # Two whitespaces are required to not find the AVX MFLOP/s by accident
        ["AVX DP MFLOP/s",              "Sum"],
        ["AVX512 DP MFLOP/s",           "Sum"],
        ["Memory bandwidth [MBytes/s]", "Sum"],
        ["Memory data volume [GBytes]", "Sum"],
        ["Local bandwidth [MByte/s]",   "Avg"],
        ["Local data volume [GByte]",   "Avg"],
        ["Remote bandwidth [MByte/s]",  "Avg"],
        ["Remote data volume [GByte]",  "Avg"],
        ["Total bandwidth [MByte/s]",   "Avg"],
        ["Total data volume [GByte]",   "Avg"],
        ["L3 bandwidth [MBytes/s]",     "Sum"],
        ["L3 data volume [GBytes]",     "Sum"],
        ["L3 request rate",             "Avg"],
        ["L3 miss rate",                "Avg"],
        ["L3 miss ratio",               "Avg"],
        ["L2 bandwidth [MBytes/s]",     "Sum"],
        ["L2 data volume [GBytes]",     "Sum"],
        ["L2 request rate",             "Avg"],
        ["L2 miss rate",                "Avg"],
        ["L2 miss ratio",               "Avg"],
        ["Branch misprediction rate",   "Avg"],
        ["Temperature [C]",             "Avg"],
        ["Energy [J]",                  "Sum"],
        ["Energy DRAM [J]",             "Max"]
]


counters = [
            ["FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE", "Sum"],
            ["FP_ARITH_INST_RETIRED_SCALAR_DOUBLE",      "Sum"],
            ["FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE", "Sum"]
           ]


def convertToFloat(val):
    '''
    Check if a val is a float and return float(val) or -1.0 if not
    Useful for likwid metrics which often use 'nil'
    '''
    try:
        return float(val)
    except:
        return -1.0

def column(matrix, i):
    return [row[i] for row in matrix]

def removeEmptyColumns(table,header):
    emptyColumns = []
    for col in range(0,len(header)):
        if all(item.strip()=="-1.0" for item in column(table,col)):
            emptyColumns.append(col)
    
    filteredTable  = []
    for row in table:
        filteredRow = []
        for col in range(0,len(header)):
            if col not in emptyColumns:
               filteredRow.append(row[col])
        filteredTable.append(filteredRow)

    filteredHeader = []
    for col in range(0,len(header)):
        if col not in emptyColumns:
            filteredHeader.append(header[col])
   
    return filteredTable,filteredHeader 

def removeInvariantColumns(table,header):
    '''
    Remove all columns containing the same value in every row
    of the given table.
    '''  
    invariantColumns        = collections.OrderedDict()
    invariantColumnsIndices = []

    blackList = ["run","run_time_steps"]

    for col in range(0,len(header)):
        current = column(table,col) 
        if header[col].strip() not in blackList\
            and all(item.strip()==current[0].strip() for item in current):
            invariantColumnsIndices.append(col)
            invariantColumns[header[col]]=current[0]
    
    filteredTable  = []
    for row in table:
        filteredRow = []
        for col in range(0,len(header)):
            if col not in invariantColumnsIndices:
               filteredRow.append(row[col])
        filteredTable.append(filteredRow)

    filteredHeader = []
    for col in range(0,len(header)):
        if col not in invariantColumnsIndices:
            filteredHeader.append(header[col])
   
    return filteredTable,filteredHeader,invariantColumns


def parseResultFile(filePath):
    '''
    Reads a single sweep job output file and parses the user time spent within each adapter.

    Args:
       filePath (str):
          Path to the sweep job output file.

    Returns:
       A dict holding for each of the found adapters a nested dict that holds the following key-value pairs:
          * 'n'       : (int)    Number of times this adapter was used.
          * 'cputime' : (float) Total CPU time spent within the adapter.
          * 'realtime': (float) Total user time spent within the adapter.

       The dict further holds the following dictionaries:
          * 'environment':(dictionary(str,str)) Total user time spent within the adapter.
          * 'parameters' :(dictionary(str,str)) Total user time spent within the adapter.
    '''
    environmentDict = {}
    parameterDict   = {}
   
    print(filePath)
 
    stats = {}
    stats["run_time_steps"]  = 0
    stats["inner_cells_min"] = 10**20
    stats["inner_cells_max"] = 0
    stats["inner_cells_avg"] = 0.0
    stats["unrefined_inner_cells_min"] = 10**20
    stats["unrefined_inner_cells_max"] = 0
    stats["unrefined_inner_cells_avg"] = 0.0
    stats["communication_time_total"] = 0.0
    stats["communication_occurences"] = 0
    stats["communication_time_avg"]   = 0.0
    stats["mesh_refinements"]     = 0.0   
    stats["local_recomputations"] = 0.0   
    stats["predictor_reruns"]     = 0.0   

    occurrences = 0   
 
    adapters      = {}
    cputimeIndex  = 3
    realtimeIndex = 5
    
    isPassedGridSetup = False
    
    # Simulating scanf functionality
    # https://docs.python.org/3/library/re.html#simulating-scanf
    commTimeRegex = re.compile("time=(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?")

    try:
        fileHandle=codecs.open(filePath,'r','UTF_8')
        for line in fileHandle:
            
            if line.startswith("sweep/environment"):
                value = line.replace("sweep/environment=","")
                environmentDict=json.loads(value)
            if line.startswith("sweep/parameters"):
                value = line.replace("sweep/parameters=","")
                parameterDict=json.loads(value)
            # inner cells/inner unrefined cells=2.13686e+06/1.81296e+06
            if "inner cells" in line:
                occurrences += 1
                m = re.search("inner cells\/inner unrefined cells=(([0-9]|\.|\+|e|E)+)\/(([0-9]|\.|\+|e|E)+)",line)
                innerCells          = float(m.group(1))
                unrefinedInnerCells = float(m.group(3))
                stats["inner_cells_min"]  = min( stats["inner_cells_min"], innerCells )
                stats["inner_cells_max"]  = max( stats["inner_cells_max"], innerCells )
                stats["inner_cells_avg"] += innerCells
                stats["unrefined_inner_cells_min"]  = min( stats["unrefined_inner_cells_min"], unrefinedInnerCells )
                stats["unrefined_inner_cells_max"]  = max( stats["unrefined_inner_cells_max"], unrefinedInnerCells )
                stats["unrefined_inner_cells_avg"] += unrefinedInnerCells
            # 154.448      [i22r02c02s11],rank:0 info         exahype::runners::Runner::startNewTimeStep(...)         step 20	t_min          =0.0015145
            m = re.search("step(\s*)([0-9]+)(\s*)t_min",line)
            if m:
                isPassedGridSetup = True
                stats["run_time_steps"] = max(stats["run_time_steps"],float(m.group(2)))
                
            # 51.4252      [cn6073.hpc.dur.ac.uk],rank:11 info         peano::performanceanalysis::DefaultAnalyser::endToPrepareAsynchronousHeapDataExchange() time=8.35e-07, cpu time=0
            if isPassedGridSetup and "endToPrepareAsynchronousHeapDataExchange()" in line:
                stats["communication_occurences"] += 1
                stats["communication_time_total"] += float(commTimeRegex.search(line).group(0).split("=")[1])
           
            # ex: nuber of mesh refinements = 0 
            if "number of mesh refinements" in line:
                stats["mesh_refinements"]     = line.split("=")[1].strip()
            if "number of local recomputations" in line:
                stats["local_recomputations"] = line.split("=")[1].strip()
            if "number of predictor reruns" in line:
                stats["predictor_reruns"]     = line.split("=")[1].strip()
 
            anchor = '|'
            header = '||'
            if anchor in line and header not in line:
                print (line)
                segments = line.split('|')
                adapter = segments[1].strip();
                adapters[adapter]                   = {}
                adapters[adapter]['iterations']     = segments[2].strip()
                adapters[adapter]['total_cputime']  = segments[cputimeIndex ].strip()
                adapters[adapter]['total_realtime'] = segments[realtimeIndex].strip()
    except IOError as err:
        print ("ERROR: could not parse adapter times for file "+filePath+"! Reason: "+str(err))
    except json.decoder.JSONDecodeError as err:
        print ("ERROR: could not parse adapter times for file "+filePath+"! Reason: "+str(err))

    if occurrences>0:
        stats["inner_cells_avg"]           = stats["inner_cells_avg"] / occurrences
        stats["unrefined_inner_cells_avg"] = stats["unrefined_inner_cells_avg"] / occurrences
    
    if stats["communication_occurences"] > 0:
        stats["communication_time_avg"] =  stats["communication_time_total"] / stats["communication_occurences"]

    return environmentDict,parameterDict,adapters,stats

def getAdapterTimesSortingKey(row):
    keyTuple = ()
    keys = row[:lastKeyColumn]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseAdapterTimes(resultsFolderPath,projectName,compressTable):
    """
    Loop over all ".out" files in the results section and create a table.
    """
    global lastKeyColumn
    lastKeyColumn=-1

    tablePath = resultsFolderPath+"/"+projectName+'.csv'
    try:
        with open(tablePath, "w") as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",", quotechar="\"")
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out")]

            unfinishedRuns = []
            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-t1-c1-b1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-N([0-9]+)-t([0-9]+)-c([0-9]+)-b([0-9]+)-r([0-9]+).out$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                ranks               = match.group(4)
                nodes               = match.group(5)
                ranksPerNode        = match.group(6)
                cores               = match.group(7)
                consumerTasks       = match.group(8)
                run                 = match.group(9)

                # derived
                totalCores = str(int(ranks)*int(cores)) # total used CPU cores
               
                environmentDict,parameterDict,adapters,stats = parseResultFile(resultsFolderPath + "/" + fileName)
                if len(environmentDict):
                    # write header
                    if firstFile:
                        header = []
                        header += sorted(environmentDict)
                        for parameter in knownParameters:
                            header.append(parameter)
                        for parameter in sorted(parameterDict):
                            if parameter not in knownParameters:
                                header.append(parameter)
                        header.append("totalCores")
                        header.append("ranks")
                        header.append("nodes")
                        header.append("ranksPerNode")
                        header.append("cores")
                        header.append("consumerTasks")
                        header.append("run")
                        header.append("adapter")
                        #
                        lastKeyColumn = len(header)
                        #
                        header.append("iterations")
                        header.append("total_cputime")
                        header.append("total_realtime")
                        header.append("normalised_cputime")
                        header.append("normalised_realtime")
                        header.append("unrefined_inner_cells_min")
                        header.append("unrefined_inner_cells_max")
                        header.append("unrefined_inner_cells_avg")
                        header.append("inner_cells_min")
                        header.append("inner_cells_max")
                        header.append("inner_cells_avg")
                        header.append("run_time_steps")
                        header.append("communication_time_total")
                        header.append("mesh_refinements")   
                        header.append("local_recomputations")
                        header.append("predictor_reruns")
                        header.append("file")
                        csvwriter.writerow(header)
                        firstFile=False

                if len(adapters): 
                    print(resultsFolderPath+"/"+fileName)

                    # write rows
                    for adapter in adapters:
                        row=[]
                        for key in sorted(environmentDict):
                            row.append(environmentDict[key])
                        for key in knownParameters:
                            row.append(parameterDict[key])
                        for key in sorted(parameterDict):
                            if key not in knownParameters:
                                row.append(parameterDict[key])
                        row.append(totalCores)
                        row.append(ranks)
                        row.append(nodes)
                        row.append(ranksPerNode)
                        row.append(cores)
                        row.append(consumerTasks)
                        row.append(run)
                        row.append(adapter)
                        row.append(adapters[adapter]["iterations"])
                        row.append(adapters[adapter]["total_cputime"])
                        row.append(adapters[adapter]["total_realtime"])
                        
                        base = -1.0
                        if "patchSize" in parameterDict:
                            base = float(parameterDict["patchSize"]) 
                        elif "order" in parameterDict:
                            base = float(parameterDict["order"]) + 1.0 
                        else:
                            print("WARNING: Found neither 'order' nor 'patchSize' key in parameter list="+str(parameterDict),file=sys.stderr)

                        print (base, parameterDict["dimension"], stats["unrefined_inner_cells_max"])
                        normalisationPerCells =  base**int(parameterDict["dimension"]) * float(stats["unrefined_inner_cells_max"])
                        if normalisationPerCells > 0:
                          row.append(str( float(adapters[adapter]["total_cputime"])  / normalisationPerCells ))
                          row.append(str( float(adapters[adapter]["total_realtime"]) / normalisationPerCells ))
                        else:
                          print("WARNING: Cannot compute normalised times as normalisation is not greater than zero. Probably, Peano's grid statistics were filtered out by a log filter.",file=sys.stderr)
                          row.append("nan")
                          row.append("nan")
                        row.append(str( int(stats["unrefined_inner_cells_min"]) ))
                        row.append(str( int(stats["unrefined_inner_cells_max"]) ))
                        row.append(str( stats["unrefined_inner_cells_avg"] ))
                        row.append(str( int(stats["inner_cells_min"]) ))
                        row.append(str( int(stats["inner_cells_max"]) ))
                        row.append(str( stats["inner_cells_avg"] ))
                        row.append(str( stats["run_time_steps"] ))
                        row.append(str( stats["communication_time_total"] ))
                        row.append(str( stats["mesh_refinements"]    ))   
                        row.append(str( stats["local_recomputations"]))   
                        row.append(str( stats["predictor_reruns"]    ))   
                        row.append(fileName)
                        csvwriter.writerow(row)
                else:
                    try: 
                        row=[]
                        for key in sorted(environmentDict):
                            row.append(environmentDict[key])
                        for key in knownParameters:
                            row.append(parameterDict[key])
                        for key in sorted(parameterDict):
                            if key not in knownParameters:
                                row.append(parameterDict[key])
                        row.append(totalCores)
                        row.append(ranks)
                        row.append(nodes)
                        row.append(ranksPerNode)
                        row.append(cores)
                        row.append(consumerTasks)
                        row.append(run)
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append("missing")
                        row.append(fileName)
                        csvwriter.writerow(row)
                    except:
                        pass
                    unfinishedRuns.append(resultsFolderPath+"/"+fileName)

        if len(unfinishedRuns):
            print("output files of unfinished runs:")
            for job in unfinishedRuns:
                print(job)

        foundFiles = not firstFile
        if foundFiles:
            # reopen the file and sort it
            tableFile   = open(tablePath, 'r')
            header      = next(tableFile)
            header      = header.strip().split(",")
            reader      = csv.reader(tableFile,delimiter=",",quotechar="\"")

            sortedData = sorted(reader,key=getAdapterTimesSortingKey)
            tableFile.close()
          
            if compressTable:
                print("") 
                sortedData,header,invariantColumns = removeInvariantColumns(sortedData,header)
                print("stripped table from the following columns as their value is the same in every row (<column header>: <common value>):")
                for column in invariantColumns:
                    print(column+": "+ invariantColumns[column])
                print("") 

            with open(tablePath, 'w') as sortedTableFile:
                writer = csv.writer(sortedTableFile, delimiter=",",quotechar="\"")
                writer.writerow(header)
                writer.writerows(sortedData)
            print("created table:")
            print(tablePath)

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: " + str(err))

def linesAreIdenticalUpToIndex(line,previousLine,index):
    result=True
    for column in range(0,index):
        result = result and line[column]==previousLine[column]
    return result

def parseSummedTimes(resultsFolderPath,projectName,timePerTimeStep=False):
    """
    Read in the sorted adapter times table and
    extract the time per time step for the fused
    and nonfused scheme.
    """
    adaptersTablePath = resultsFolderPath+"/"+projectName+'.csv'
    outputTablePath   = resultsFolderPath+"/"+projectName+'-total-times.csv'
    if timePerTimeStep:
        outputTablePath   = resultsFolderPath+"/"+projectName+'-timestep-times.csv'
    try:
        print ("reading table "+adaptersTablePath+"")
        adaptersTableFile  = open(adaptersTablePath, 'r')
        header             = next(adaptersTableFile).strip().split(',')
        csvreader          = csv.reader(adaptersTableFile,delimiter=",",quotechar="\"")

        runColumn                = header.index("run")
        adapterColumn            = header.index("adapter")
        iterationsColumn         = header.index("iterations")
        cpuTimeColumn            = header.index("total_cputime")
        userTimeColumn           = header.index("total_realtime")
        normalisedCPUTimeColumn  = header.index("normalised_cputime")
        normalisedUserTimeColumn = header.index("normalised_realtime")
        runTimeStepsColumn       = header.index("run_time_steps")
        
        if runColumn >= adapterColumn:
            print ("ERROR: order of columns not suitable. Column 'run' must come before column 'adapter'!")
        
        def appendMoments(row,measurements):
            row.append(str(min(measurements)))
            row.append(str(max(measurements)))
            row.append(str(statistics.mean(measurements)))
            if len(measurements)>1:
                row.append(str(statistics.stdev(measurements)))
            else:
                row.append("0.0")
        
        fused = True
        
        def sumUpTimes(line):
            nonlocal fused
        
            fusedAdapters        = ["FusedTimeStep"]
            nonfusedAdapters     = ["Correction","MergeNeighbours","Prediction","UpdateAndReduce"]
            
            adapter = line[adapterColumn]
            if adapter in fusedAdapters:
              fused = True
            elif adapter in ["Correction", "MergeNeighbours"]:
              fused = False
            
            if timePerTimeStep and (fused and adapter in fusedAdapters) or (not fused and adapter in nonfusedAdapters):
                summedCPUTimes[-1]            += float(line[cpuTimeColumn]) / float(line[runTimeStepsColumn])
                summedUserTimes[-1]           += float(line[userTimeColumn]) / float(line[runTimeStepsColumn])
                summedNormalisedCPUTimes[-1]  += float(line[normalisedCPUTimeColumn]) / float(line[runTimeStepsColumn])
                summedNormalisedUserTimes[-1] += float(line[normalisedUserTimeColumn]) / float(line[runTimeStepsColumn])
            elif not timePerTimeStep:
                summedCPUTimes[-1]            += float(line[cpuTimeColumn])
                summedUserTimes[-1]           += float(line[userTimeColumn])
                summedNormalisedCPUTimes[-1]  += float(line[normalisedCPUTimeColumn])
                summedNormalisedUserTimes[-1] += float(line[normalisedUserTimeColumn])
       
        # open file 
        with open(outputTablePath, 'w') as timeStepTimesTableFile:
            csvwriter = csv.writer(timeStepTimesTableFile, delimiter=",",quotechar="\"")
            # write header
            row = header[0:runColumn]
            row.append("runs")
            row.append("cputime_min")
            row.append("cputime_max")
            row.append("cputime_mean")
            row.append("cputime_stdev")
            row.append("realtime_min")
            row.append("realtime_max")
            row.append("realtime_mean")
            row.append("realtime_stdev")
            row.append("normalised_cputime_min")
            row.append("normalised_cputime_max")
            row.append("normalised_cputime_mean")
            row.append("normalised_cputime_stdev")
            row.append("normalised_realtime_min")
            row.append("normalised_realtime_max")
            row.append("normalised_realtime_mean")
            row.append("normalised_realtime_stdev")
            csvwriter.writerow(row)
            
            # init
            summedCPUTimes            = [0.0]
            summedUserTimes           = [0.0]
            summedNormalisedCPUTimes  = [0.0]
            summedNormalisedUserTimes = [0.0]
            
            previousLine      = None
            # write intermediate rows
            for line in csvreader:
                if previousLine==None:
                    previousLine=line
                
                # new adapter 
                if linesAreIdenticalUpToIndex(line,previousLine,adapterColumn):
                    adapter = line[adapterColumn]
                    #print("new adapter: "+adapter)
                    if adapter!="missing":
                        sumUpTimes(line)
                # new run  
                elif linesAreIdenticalUpToIndex(line,previousLine,runColumn): 
                    adapter = line[adapterColumn]
                    #print("new run: "+adapter)
                    if adapter!="missing":
                        summedCPUTimes.append(0.0)
                        summedUserTimes.append(0.0)
                        summedNormalisedCPUTimes.append(0.0)
                        summedNormalisedUserTimes.append(0.0)
 
                        sumUpTimes(line)
                # new experiment
                else:
                    row = previousLine[0:runColumn]
                    row.append(str(len(summedCPUTimes)))

                    if len(summedCPUTimes): 
                        appendMoments(row,summedCPUTimes)
                        appendMoments(row,summedUserTimes)
                        appendMoments(row,summedNormalisedCPUTimes)
                        appendMoments(row,summedNormalisedUserTimes)
                    else:
                        summedCPUTimes            = [float("nan")]
                        summedUserTimes           = [float("nan")]
                        summedNormalisedCPUTimes  = [float("nan")]
                        summedNormalisedUserTimes = [float("nan")]
                        appendMoments(row,summedCPUTimes)
                        appendMoments(row,summedUserTimes)
                        appendMoments(row,summedNormalisedCPUTimes)
                        appendMoments(row,summedNormalisedUserTimes)
                    
                    csvwriter.writerow(row)
                    # reset 
                    summedCPUTimes            = [0.0]
                    summedUserTimes           = [0.0]
                    summedNormalisedCPUTimes  = [0.0]
                    summedNormalisedUserTimes = [0.0]
                    
                    adapter = line[adapterColumn]
                    # print("new experiment: "+adapter)
                    if adapter!="missing":
                        summedCPUTimes            = [0.0]
                        summedUserTimes           = [0.0]
                        summedNormalisedCPUTimes  = [0.0]
                        summedNormalisedUserTimes = [0.0]
                   
                        sumUpTimes(line)
                    else:
                        summedCPUTimes            = []
                        summedUserTimes           = []
                        summedNormalisedCPUTimes  = []
                        summedNormalisedUserTimes = []
                
                previousLine  = line
            
            # write last row (copy and paste)
            row = previousLine[0:runColumn]
            row.append(str(len(summedCPUTimes)))
            if len(summedCPUTimes): 
                appendMoments(row,summedCPUTimes)
                appendMoments(row,summedUserTimes)
                appendMoments(row,summedNormalisedCPUTimes)
                appendMoments(row,summedNormalisedUserTimes)
            else:
                summedCPUTimes            = [float("nan")]
                summedUserTimes           = [float("nan")]
                summedNormalisedCPUTimes  = [float("nan")]
                summedNormalisedUserTimes = [float("nan")]
                appendMoments(row,summedCPUTimes)
                appendMoments(row,summedUserTimes)
                appendMoments(row,summedNormalisedCPUTimes)
                appendMoments(row,summedNormalisedUserTimes)     
            csvwriter.writerow(row)
            
            print("created table:")
            print(outputTablePath)
    except IOError as err:
        print ("ERROR: could not read file "+adaptersTablePath+". Error message: "<<str(err))


def parseLikwidMetrics(filePath,metrics,counters,singlecore=False):
    """
    Reads a single Peano output file and parses likwid performance metrics.

    Args:
       filePath (str):
          Path to the Peano output file.
       metrics (str[][]):
          A list of metrics the we want to read out.
       counters (str[][]):
          A list of counters the we want to read out.
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

    environmentDict = {}
    parameterDict   = {}

    result  = {}
    for metric in metrics:
        result[metric[0]] =  {}
        result[metric[0]][metric[1]] = -1.0
    for counter in counters:
        result[counter[0]] =  {}
        result[counter[0]][counter[1]] = -1.0

    try:
        fileHandle=open(filePath,encoding="utf-8")

        for line in fileHandle:
            if line.startswith("sweep/environment"):
                value = line.split('=')[-1]
                environmentDict=json.loads(value)
            if line.startswith("sweep/parameters"):
                value = line.split('=')[-1]
                parameterDict=json.loads(value)

            for metric in metrics:
                if metric[0]+" STAT" in line:
                    segments = line.split('|')
                    #   |  Runtime (RDTSC) [s] STAT |   27.4632  |   1.1443  |   1.1443  |   1.1443  |
                    values = {}
                    values["Sum"] = convertToFloat(segments[2].strip());
                    values["Min"] = convertToFloat(segments[3].strip());
                    values["Max"] = convertToFloat(segments[4].strip());
                    values["Avg"] = convertToFloat(segments[5].strip());
                    result[metric[0]][metric[1]]=values[metric[1]]
                elif metric[0] in line:
                    segments = line.split('|')

                    #    |     Runtime (RDTSC) [s]    |    6.5219    |
                    value  = convertToFloat(segments[3].strip());
                    values = {}
                    values["Sum"] = value
                    values["Min"] = value
                    values["Max"] = value
                    values["Avg"] = value
                    result[metric[0]][metric[1]]=values[metric[1]]

            for counter in counters:
                if counter[0]+" STAT" in line:
                    segments = line.split('|')
                    #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  623010105225  | ...
                    values = {}
                    values["Sum"] = convertToFloat(segments[3].strip())
                    values["Min"] = convertToFloat(segments[4].strip())
                    values["Max"] = convertToFloat(segments[5].strip())
                    values["Avg"] = convertToFloat(segments[6].strip())
                    result[counter[0]][counter[1]]=values[counter[1]]
                elif counter[0] in line:
                    segments = line.split('|')
                    #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  623010105225  | ...
                    value  = convertToFloat(segments[3].strip());
                    values = {}
                    values["Sum"] = value
                    values["Min"] = value
                    values["Max"] = value
                    values["Avg"] = value
                    result[counter[0]][counter[1]]=values[metric[1]]

    except IOError as err:
        print ("ERROR: could not parse likwid metrics for file "+filePath+"! Reason: "+str(err))
    return environmentDict,parameterDict,result

def getLikwidMetricsSortingKey(row):
    keyTuple = ()
    keys = row[:-(len(metrics)+len(counters))]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseMetrics(resultsFolderPath,projectName,compressTable):
    """
    Loop over all ".out.likwid" files in the results section and create a table.
    """
    
    tablePath         = resultsFolderPath+"/"+projectName+'-likwid.csv'
    try:
        with open(tablePath, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",",quotechar="\"")
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out.likwid")]
            
            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-N1-t1-c1-b1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-N([0-9]+)-t([0-9]+)-c([0-9]+)-b([0-9]+)-r([0-9]+).out.likwid$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                ranks               = match.group(4)
                nodes               = match.group(5)
                ranksPerNode        = match.group(6)
                cores               = match.group(7)
                consumerTasks       = match.group(8)
                run                 = match.group(9)

                environmentDict,parameterDict,measurements = parseLikwidMetrics(resultsFolderPath + "/" + fileName, metrics, counters, cores.startswith("1:"))

                # TODO(Dominic): workaround. parameters
                if len(environmentDict) is 0:
                   environmentDict,parameterDict,adapters,stats = parseResultFile(resultsFolderPath + "/" + fileName.replace(".likwid",""))

                if len(measurements):
                    # write header
                    if firstFile:
                        header = []
                        header += sorted(environmentDict)
                        for parameter in knownParameters:
                            header.append(parameter)
                        for parameter in sorted(parameterDict):
                            if parameter not in knownParameters:
                                header.append(parameter)
                        header.append("ranks")
                        header.append("nodes")
                        header.append("ranksPerNode")
                        header.append("cores")
                        header.append("consumerTasks")
                        header.append("run")
                        for key in sorted(measurements):
                            for subkey in measurements[key]:
                                header.append(key+" ("+subkey+")")
                        header.append("file")
                        csvwriter.writerow(header)
                        firstFile=False
                    print(resultsFolderPath+"/"+fileName)

                    # write row
                    row=[]
                    for key in sorted(environmentDict):
                        row.append(environmentDict[key])
                    for key in knownParameters:
                        row.append(parameterDict[key])
                    for key in sorted(parameterDict):
                        if key not in knownParameters:
                            row.append(parameterDict[key])
                    row.append(ranks)
                    row.append(nodes)
                    row.append(ranksPerNode)
                    row.append(cores)
                    row.append(consumerTasks)
                    row.append(run)
                    for key in sorted(measurements):
                        for subkey in measurements[key]:
                            row.append(measurements[key][subkey])
                    row.append(fileName)
                    csvwriter.writerow(row)

        success = not firstFile
        if success:
          # reopen the table and sort it
          tableFile   = open(tablePath, 'r')
          header      = next(tableFile)
          header      = header.strip().split(",")
          reader      = csv.reader(tableFile,delimiter=",",quotechar="\"")

          sortedData = sorted(reader,key=getLikwidMetricsSortingKey)
          tableFile.close()

          if compressTable:
             print("") 
             sortedData,header = removeEmptyColumns(sortedData,header)
             print("stripped table from columns containing value \"-1.0\" in every row.") 
             sortedData,header,invariantColumns = removeInvariantColumns(sortedData,header)
             print("stripped table from the following columns as their value is the same in every row (<column header>: <common value>):")
             for column in invariantColumns:
                 print(column+": "+ invariantColumns[column])
             print("") 

          with open(tablePath, 'w') as sortedTableFile:
              writer = csv.writer(sortedTableFile, delimiter=",",quotechar="\"")
              writer.writerow(header)
              writer.writerows(sortedData)
          print("created table:")
          print(tablePath)

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: "<<str(err))

def parseJobStatistics(resultsFolderPath,projectName,compressTable):
    """
    Loop over all ".out" files in the results section and parses the
    job statistics.
    """
    global lastKeyColumn
    lastKeyColumn=-1
    
    tablePath = resultsFolderPath+"/"+projectName+'-job-statistics.csv'
    try:
        with open(tablePath, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",",quotechar="\"")
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out")]
            
            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-N1-t1-c1-b1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-N([0-9]+)-t([0-9]+)-c([0-9]+)-b([0-9]+)-r([0-9]+).out$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                ranks               = match.group(4)
                nodes               = match.group(5)
                ranksPerNode        = match.group(6)
                cores               = match.group(7)
                consumerTasks       = match.group(8)
                run                 = match.group(9)

                environmentDict,parameterDict,statsNoOfBackgroundTasks,\
                statsGrabbedBackgroundTasks,statsRunningConsumers =\
                    parseJobStatisticsFromResultsFile(resultsFolderPath + "/" + fileName)

                if environmentDict!=None:
                    # write header
                    if firstFile:
                        header = []
                        header += sorted(environmentDict)
                        for parameter in knownParameters:
                            header.append(parameter)
                        for parameter in sorted(parameterDict):
                            if parameter not in knownParameters:
                                header.append(parameter)
                        header.append("ranks")
                        header.append("nodes")
                        header.append("ranksPerNode")
                        header.append("cores")
                        header.append("consumerTasks")
                        header.append("run")
                        #
                        lastKeyColumn = len(header)
                        #
                        for bucket in sorted(statsNoOfBackgroundTasks):
                            header.append("numJobs{}".format(bucket))
                        for bucket in sorted(statsGrabbedBackgroundTasks):
                            header.append("numGrabbed{}".format(bucket))
                        for bucket in sorted(statsRunningConsumers):
                            header.append("numConsumers{}".format(bucket))
                        header.append("file")
                        csvwriter.writerow(header)
                        firstFile=False
                    print(resultsFolderPath+"/"+fileName)

                    # write row
                    row=[]
                    for key in sorted(environmentDict):
                        row.append(environmentDict[key])
                    for key in knownParameters:
                        row.append(parameterDict[key])
                    for key in sorted(parameterDict):
                        if key not in knownParameters:
                            row.append(parameterDict[key])
                    row.append(ranks)
                    row.append(nodes)
                    row.append(ranksPerNode)
                    row.append(cores)
                    row.append(consumerTasks)
                    row.append(run)
                    for bucket in sorted(statsNoOfBackgroundTasks):
                        row.append(statsNoOfBackgroundTasks[bucket])
                    for bucket in sorted(statsGrabbedBackgroundTasks):
                        row.append(statsGrabbedBackgroundTasks[bucket])
                    for bucket in sorted(statsRunningConsumers):
                        row.append(statsRunningConsumers[bucket])
                    row.append(fileName)
                    csvwriter.writerow(row)

        success = not firstFile
        if success:
          # reopen the table and sort it
          tableFile   = open(tablePath, 'r')
          header      = next(tableFile)
          header      = header.strip().split(",")
          reader      = csv.reader(tableFile,delimiter=",",quotechar="\"")

          sortedData = sorted(reader,key=getLikwidMetricsSortingKey)
          tableFile.close()

          if compressTable:
             print("") 
             sortedData,header = removeEmptyColumns(sortedData,header)
             print("stripped table from columns containing value \"-1.0\" in every row.") 
             sortedData,header,invariantColumns = removeInvariantColumns(sortedData,header)
             print("stripped table from the following columns as their value is the same in every row (<column header>: <common value>):")
             for column in invariantColumns:
                 print(column+": "+ invariantColumns[column])
             print("") 

          with open(tablePath, 'w') as sortedTableFile:
              writer = csv.writer(sortedTableFile, delimiter=",",quotechar="\"")
              writer.writerow(header)
              writer.writerows(sortedData)
          print("created table:")
          print(tablePath)

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: "<<str(err))
    return

def parseJobStatisticsFromResultsFile(resultsFile):
    statsNoOfBackgroundTasks = {}
    statsGrabbedBackgroundTasks = {}
    statsRunningConsumers = {}
    environmentDict = None
    parameterDict = None   
 
    # bins: < x< 2^0, x< 2^1, x < 2^2
    for bucket in range(0,29): # ~ 270 milloon tasks is largest bin 2^28
        statsNoOfBackgroundTasks[bucket]    = 0
        statsGrabbedBackgroundTasks[bucket] = 0
    for bucket in range(0,11): # 1024 consumers can run on the same time
        statsRunningConsumers[bucket] = 0
    #  10.4573      info         no of background tasks[1]=19029 
    base = 2
    with open(resultsFile, 'r') as f:
        for l in f.readlines():
            if l.startswith("sweep/environment"):
                value = l.replace("sweep/environment=","")
                environmentDict=json.loads(value)
            if l.startswith("sweep/parameters"):
                value = l.replace("sweep/parameters=","")
                parameterDict=json.loads(value)
            if "no of background tasks available per consumer run" in l:
                ex = re.compile(r"no of background tasks available per consumer run\[(\d+)\]=(\d+)?")
                match = ex.search(l)
                if match:
                    num_tasks = int(match.group(1))
                    occurences = int(match.group(2))
                    if num_tasks == 0:
                        statsNoOfBackgroundTasks[0] += occurences
                    else:
                        bucket = int(math.log(num_tasks, base)) + 1
                        statsNoOfBackgroundTasks[bucket] += occurences
            if "no of background tasks processed per consumer run" in l:
                ex = re.compile(r"no of background tasks processed per consumer run\[(\d+)\]=(\d+)?")
                match = ex.search(l)
                if match:
                    num_tasks = int(match.group(1))
                    occurences = int(match.group(2))
                    if num_tasks == 0:
                        statsGrabbedBackgroundTasks[0] += occurences
                    else:
                        bucket = int(math.log(num_tasks, base)) + 1
                        statsGrabbedBackgroundTasks[bucket] += occurences
            if "no of running consumers" in l:
                ex = re.compile(r"no of running consumers\[(\d+)\]=(\d+)?")
                match = ex.search(l)
                if match:
                    num_tasks  = int(match.group(1))
                    occurences = int(match.group(2))
                    if num_tasks == 0:
                        statsRunningConsumers[0] += occurences
                    else:
                        bucket = int(math.log(num_tasks, base)) + 1
                        statsRunningConsumers[bucket] += occurences
    return environmentDict,parameterDict,statsNoOfBackgroundTasks,statsGrabbedBackgroundTasks,statsRunningConsumers


def parseArgs():
    parser = argparse.ArgumentParser(
        description="A collection of parsers for analysing ExaHyPE application output.",
            epilog="End of help message.",
    )
    parser.add_argument("--compress", dest="compress", action="store_true",help="Remove columns where the same value is found in every row.")
    parser.set_defaults(compress=False)

    parser.add_argument("--prefix",nargs="?",default=None, help="Specify the prefix of the sweep/ExaHyPE output files.")
    parser.set_defaults(prefix=None)
    
    # subprograms
    parser.add_argument("--parseAllMetrics", dest="parseAllMetrics", action="store_true",help="Parse all output files in the specified directory with the specified prefix.")
    parser.set_defaults(parseAllMetrics=False)
    parser.add_argument("--parseAllAdapters", dest="parseAllAdapters", action="store_true",help="Parse all output files in the specified directory with the specified prefix.")
    parser.set_defaults(parseAllAdapters=False)
    parser.add_argument("--parseAllTimeStepTimes", dest="parseAllTimeStepTimes", action="store_true",help="Parse all output files in the specified directory with the specified prefix.")
    parser.set_defaults(parseAllTimeStepTimes=False)
    
    parser.add_argument("file",
        type=str,help="The directory or CSV file to work with.")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parseArgs();

    if args.parseAllMetrics and args.prefix!=None and os.path.isdir(args.file):
        parseMetrics(args.file,args.prefix,args.compress)
    if args.parseAllAdapters and args.prefix!=None and os.path.isdir(args.file):
        parseAdapterTimes(args.file,args.prefix,args.compress)
    if args.parseAllTimeStepTimes and args.prefix!=None and os.path.isdir(args.file):
        parseSummedTimes(args.file,args.prefix,True)
