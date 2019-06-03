#!/usr/bin/env python3
"""
.. module:: sweep_analysis
  :platform: Unix, Windows, Mac
  :synopsis: parse L1, L2, and LInf error from output files
   and creating

.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>,

:synopsis: Generate benchmark suites for ExaHyPE.
"""
import csv
import json
import sys,os
import re

import statistics

import argparse

knownParameters = ["dimension","architecture"]

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

def getSortingKey(row):
    keyTuple = ()
    keys = row[:lastKeyColumn]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseReducedValues(resultsFolderPath,projectName,compressTable):
    """
    Loop over all ".out" files in the specified folder and parses
    data that has the following form: '/sweep/reduce/<key>=<value>'.

    For every found key, the subroutine will compute min, max, mean 
    and standard deviation for all values in the output file.
    """
    global lastKeyColumn
    lastKeyColumn=-1
    
    tablePath = resultsFolderPath+"/"+projectName+'-reductions.csv'
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

                environmentDict,parameterDict,reductions = \
                    parseReducedValuesFromResultsFile(resultsFolderPath + "/" + fileName)

                if environmentDict!=None and len(reductions):
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
                        for key in sorted(reductions):
                            for suffix in ["min","max","avg","stdev"]:
                                header.append("%s_%s" % (key,suffix))
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
                    for key in sorted(reductions):
                         row.append(min(reductions[key]))
                         row.append(max(reductions[key]))
                         row.append(statistics.mean(reductions[key]))
                         row.append(statistics.stdev(reductions[key]))
                    row.append(fileName)
                    csvwriter.writerow(row)

        success = not firstFile
        if success:
          # reopen the table and sort it
          tableFile   = open(tablePath, 'r')
          header      = next(tableFile)
          header      = header.strip().split(",")
          reader      = csv.reader(tableFile,delimiter=",",quotechar="\"")

          sortedData = sorted(reader,key=getSortingKey)
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

def parseReducedValuesFromResultsFile(resultsFile): 
    reductions = {} 
    environmentDict = None
    parameterDict = None   
 
    with open(resultsFile, 'r') as f:
        for l in f.readlines():
            if l.startswith("sweep/environment"):
                value = l.replace("sweep/environment=","")
                environmentDict=json.loads(value)
            if l.startswith("sweep/parameters"):
                value = l.replace("sweep/parameters=","")
                parameterDict=json.loads(value)
            if "sweep/reduce/" in l:
                stripped=value = l.replace("sweep/reduce/","")
                tokens=stripped.split("=")
                if len(tokens)==2:
                  key   = tokens[0].split(" ")[-1].strip()
                  value = tokens[1].strip()
                  try:
                      value=float(value)
                  except:
                      value=float("nan")
                  if key in reductions:
                      reductions[key].append(value)
                  else:
                      reductions[key] = [value]

    # stddev computation requires at least two values
    for key in reductions:
        if len(reductions[key])==1:
             reductions[key].append(reductions[key][0])

    return environmentDict,parameterDict,reductions

def parseArgs():
    parser = argparse.ArgumentParser(
        description="A collection of parsers for analysing ExaHyPE application output.",
            epilog="End of help message.",
    )
    parser.add_argument("--compress", dest="compress", action="store_true",help="Remove columns where the same value is found in every row.")
    parser.set_defaults(compress=False)

    parser.add_argument("--prefix",nargs="?",default=None, help="Specify the prefix of the sweep/ExaHyPE output files.")
    parser.set_defaults(prefix=None)
    
    parser.add_argument("file",
        type=str,help="The directory or CSV file to work with.")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parseArgs();

    parseReducedValues(args.file,args.prefix,args.compress)
