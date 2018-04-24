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
import os
import re
import codecs

import statistics

knownParameters   = ["architecture", "optimisation", "dimension", "order" ]

metrics =  [
        ["  MFLOP/s",                   "Sum"],  # Two whitespaces are required to not find the AVX MFLOP/s by accident
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

counters = [
            ["FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE", "Sum"],
            ["FP_ARITH_INST_RETIRED_SCALAR_DOUBLE",      "Sum"],
            ["FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE", "Sum"]
           ]

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
          * 'usertime': (float) Total user time spent within the adapter.

       The dict further holds the following dictionaries:
          * 'environment':(dictionary(str,str)) Total user time spent within the adapter.
          * 'parameters' :(dictionary(str,str)) Total user time spent within the adapter.
    '''
    environmentDict = {}
    parameterDict   = {}
    
    stats = {}
    stats["inner_cells_min"] = 10.0**20
    stats["inner_cells_max"] = 0.0
    stats["inner_cells_avg"] = 0.0
    stats["unrefined_inner_cells_min"] = 10.0**20
    stats["unrefined_inner_cells_max"] = 0.0
    stats["unrefined_inner_cells_avg"] = 0.0
    
    adapters      = {}
    cputimeIndex  = 3
    usertimeIndex = 5

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
                m = re.search("inner cells\/inner unrefined cells=(([0-9]|\.|\+|e|E)+)\/(([0-9]|\.|\+|e|E)+)",line)
                innerCells          = float(m.group(1))
                unrefinedInnerCells = float(m.group(3))
                stats["inner_cells_min"]  = min( ["inner_cells_min"], innerCells )
                stats["inner_cells_max"]  = min( ["inner_cells_max"], innerCells )
                stats["inner_cells_avg"] += innerCells
                stats["unrefined_inner_cells_min"]  = min( ["unrefined_inner_cells_min"], innerCells )
                stats["unrefined_inner_cells_max"]  = min( ["unrefined_inner_cells_max"], innerCells )
                stats["unrefined_inner_cells_avg"] += innerCells
                
            anchor = '|'
            header = '||'
            if anchor in line and header not in line:
                segments = line.split('|')
                adapter = segments[1].strip();
                adapters[adapter]                   = {}
                adapters[adapter]['iterations']     = segments[2].strip()
                adapters[adapter]['total_cputime']  = segments[cputimeIndex ].strip()
                adapters[adapter]['total_usertime'] = segments[usertimeIndex].strip()
    except IOError as err:
        print ("ERROR: could not parse adapter times for file "+filePath+"! Reason: "+str(err))
    except json.decoder.JSONDecodeError as err:
        print ("ERROR: could not parse adapter times for file "+filePath+"! Reason: "+str(err))
    return environmentDict,parameterDict,adapters,statistics

def getAdapterTimesSortingKey(row):
    keyTuple = ()
    keys = row[:-3]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseAdapterTimes(resultsFolderPath,projectName):
    """
    Loop over all ".out" files in the results section and create a table.
    """
    tablePath = resultsFolderPath+"/"+projectName+'.csv'
    try:
        with open(tablePath, "w") as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",", quotechar="\"")
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out")]

            unfinishedRuns = []
            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-t1-c1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-N([0-9]+)-t([0-9]+)-c([0-9]+)-r([0-9]+).out$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                ranks               = match.group(4)
                nodes               = match.group(5)
                tasks               = match.group(6)
                cores               = match.group(7)
                run                 = match.group(8)
                
                environmentDict,parameterDict,adapters,stats = parseResultFile(resultsFolderPath + "/" + fileName)
                if len(adapters):
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
                        header.append("tasks")
                        header.append("cores")
                        header.append("run")
                        header.append("adapter")
                        header.append("iterations")
                        header.append("total_cputime")
                        header.append("total_usertime")
                        header.append("normalised_cputime")
                        header.append("normalised_usertime")
                        header.append("unrefined_inner_cells_min")
                        header.append("unrefined_inner_cells_max")
                        header.append("unrefined_inner_cells_avg")
                        header.append("inner_cells_min")
                        header.append("inner_cells_max")
                        header.append("inner_cells_avg")
                        header.append("file")
                        csvwriter.writerow(header)
                        firstFile=False
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
                        row.append(ranks)
                        row.append(nodes)
                        row.append(tasks)
                        row.append(cores)
                        row.append(run)
                        row.append(adapter)
                        row.append(adapters[adapter]["iterations"])
                        row.append(adapters[adapter]["total_cputime"])
                        row.append(adapters[adapter]["total_usertime"])
                        
                        normalisationPerCells =  float(adapters[adapter]["iterations"]) *( float(parameterDict["order"]) + 1 )**int(parameterDict["dimension"]) *\
                                                 stats["unrefined_inner_cells_avg"]
                        row.append(str( float(adapters[adapter]["total_cputime"])  / normalisationPerCells ))
                        row.append(str( float(adapters[adapter]["total_usertime"]) / normalisationPerCells ))
                        row.append(str( stats["unrefined_inner_cells_min"] ))
                        row.append(str( stats["unrefined_inner_cells_max"] ))
                        row.append(str( stats["unrefined_inner_cells_avg"] ))
                        row.append(str( stats["inner_cells_min"] ))
                        row.append(str( stats["inner_cells_max"] ))
                        row.append(str( stats["inner_cells_avg"] ))
                        row.append(fileName)
                        csvwriter.writerow(row)
                else:
                    unfinishedRuns.append(resultsFolderPath+"/"+fileName)
        if len(unfinishedRuns):
            print("output files of unfinished runs:")
            for job in unfinishedRuns:
                print(job)

        success = not firstFile
        if success:
            # reopen the file and sort it
            tableFile   = open(tablePath, 'r')
            header      = next(tableFile)
            header      = header.strip()
            reader      = csv.reader(tableFile,delimiter=",",quotechar="\"")

            sortedData = sorted(reader,key=getAdapterTimesSortingKey)
            tableFile.close()

            with open(tablePath, 'w') as sortedTableFile:
                writer = csv.writer(sortedTableFile, delimiter=",",quotechar="\"")
                writer.writerow(header.split(','))
                writer.writerows(sortedData)
            print("created table:")
            print(tablePath)

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: " + str(err))


def column(matrix, i):
    return [row[i] for row in matrix]

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
    fusedAdapters        = ["FusedTimeStep"]
    firstNonfusedAdapter = "BroadcastAndMergeNeighbours"
    nonfusedAdapters     = [firstNonfusedAdapter, "Prediction", "UpdateAndReduce"]

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
        userTimeColumn           = header.index("total_usertime")
        normalisedCPUTimeColumn  = header.index("normalised_cputime")
        normalisedUserTimeColumn = header.index("normalised_usertime")
        
        if runColumn >= adapterColumn:
            print ("ERROR: order of columns not suitable. Column 'run' must come before column 'adapter'!")
        
        with open(outputTablePath, 'w') as timeStepTimesTableFile:
            csvwriter = csv.writer(timeStepTimesTableFile, delimiter=",",quotechar="\"")
            # write header
            row = header[0:runColumn]
            row.append("runs")
            row.append("cputime_min")
            row.append("cputime_max")
            row.append("cputime_mean")
            row.append("cputime_stdev")
            row.append("usertime_min")
            row.append("usertime_max")
            row.append("usertime_mean")
            row.append("usertime_stdev")
            row.append("normalised_cputime_min")
            row.append("normalised_cputime_max")
            row.append("normalised_cputime_mean")
            row.append("normalised_cputime_stdev")
            row.append("normalised_usertime_min")
            row.append("normalised_usertime_max")
            row.append("normalised_usertime_mean")
            row.append("normalised_usertime_stdev")
            csvwriter.writerow(row)
            
            # init
            summedCPUTimes            = [0.0]
            summedUserTimes           = [0.0]
            summedNormalisedCPUTimes  = [0.0]
            summedNormalisedUserTimes = [0.0]
            
            def appendMoments(row,measurements):
                row.append(str(min(measurements)))
                row.append(str(max(measurements)))
                row.append(str(statistics.mean(measurements)))
                if len(measurements)>1:
                    row.append(str(statistics.stdev(measurements)))
                else:
                    row.append("0.0")
            
            previousLine      = None
            fused             = True
            # write intermediate rows
            for line in csvreader:
                if previousLine==None:
                    previousLine=line
                
                if linesAreIdenticalUpToIndex(line,previousLine,adapterColumn):
                    adapter = line[adapterColumn]
                    if adapter==firstNonfusedAdapter:
                       fused = False
                    
                    if timePerTimeStep and (fused and adapter in fusedAdapters) or (not fused and adapter in nonfusedAdapters):
                        summedCPUTimes[-1]            += float(line[cpuTimeColumn]) / float(line[iterationsColumn])
                        summedUserTimes[-1]           += float(line[userTimeColumn]) / float(line[iterationsColumn])
                        summedNormalisedCPUTimes[-1]  += float(line[normalisedCPUTimeColumn]) / float(line[iterationsColumn])
                        summedNormalisedUserTimes[-1] += float(line[normalisedUserTimeColumn]) / float(line[iterationsColumn])
                    elif not timePerTimeStep:
                        summedCPUTimes[-1]            += float(line[cpuTimeColumn])
                        summedUserTimes[-1]           += float(line[userTimeColumn])
                        summedNormalisedCPUTimes[-1]  += float(line[normalisedCPUTimeColumn])
                        summedNormalisedUserTimes[-1] += float(line[normalisedUserTimeColumn])
                elif linesAreIdenticalUpToIndex(line,previousLine,runColumn):
                    summedCPUTimes.append(0.0)
                    summedUserTimes.append(0.0)
                    summedNormalisedCPUTimes.append(0.0)
                    summedNormalisedUserTimes.append(0.0)
                else:
                    row = previousLine[0:runColumn]
                    row.append(str(len(summedCPUTimes)))
                    
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
                    
                    fused = True
                
                previousLine  = line
            
            # write last row (copy and paste)
            row = previousLine[0:runColumn]
            row.append(str(len(summedCPUTimes)))
            
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
        fileHandle=open(filePath)

        for line in fileHandle:
            if line.startswith("sweep/environment"):
                value = line.split('=')[-1]
                environmentDict=json.loads(value)
            if line.startswith("sweep/parameters"):
                value = line.split('=')[-1]
                parameterDict=json.loads(value)

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

            for counter in counters:
                if singlecore:
                    if counter[0] in line:
                        segments = line.split('|')
                        #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  623010105225  | ...
                        value  = float(segments[3].strip());
                        values = {}
                        values["Sum"] = value
                        values["Min"] = value
                        values["Max"] = value
                        values["Avg"] = value
                        result[counter[0]][counter[1]]=values[metric[1]]
                else:
                    if counter[0]+" STAT" in line:
                        segments = line.split('|')
                        #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  623010105225  | ...
                        values = {}
                        values["Sum"] = float(segments[3].strip());
                        values["Min"] = float(segments[4].strip());
                        values["Max"] = float(segments[5].strip());
                        values["Avg"] = float(segments[6].strip());
                        result[counter[0]][counter[1]]=values[counter[1]]
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

def parseMetrics(resultsFolderPath,projectName):
    """
    Loop over all ".out.likwid" files in the results section and create a table.
    """
    
    tablePath         = resultsFolderPath+"/"+projectName+'-likwid.csv'
    try:
        with open(tablePath, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=";")
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out.likwid")]
            
            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-N1-t1-c1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-N([0-9]+)-t([0-9]+)-c([0-9]+)-r([0-9]+).out.likwid$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                ranks               = match.group(4)
                nodes               = match.group(5)
                tasks               = match.group(6)
                cores               = match.group(7)
                run                 = match.group(8)

                environmentDict,parameterDict,measurements = parseLikwidMetrics(resultsFolderPath + "/" + fileName, metrics, counters, cores=="1")

                # TODO(Dominic): workaround. parameters
                if len(environmentDict) is 0:
                   environmentDict,parameterDict,adapters = parseResultFile(resultsFolderPath + "/" + fileName.replace(".likwid",""))

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
                        header.append("tasks")
                        header.append("cores")
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
                    row.append(tasks)
                    row.append(cores)
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
          header      = header.strip()
          reader      = csv.reader(tableFile,delimiter=";")

          sortedData = sorted(reader,key=getLikwidMetricsSortingKey)
          tableFile.close()

          with open(tablePath, 'w') as sortedTableFile:
              writer = csv.writer(sortedTableFile, delimiter=";")
              writer.writerow(header.split(";"))
              writer.writerows(sortedData)
          print("created table:")
          print(tablePath)

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: "<<str(err))
