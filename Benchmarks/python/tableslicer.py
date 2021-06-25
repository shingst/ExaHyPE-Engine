#!/usr/bin/env python3
def column(matrix, i):
    return [row[i] for row in matrix]


def minAndMaxInColumn(matrix, i):
     minimum = +float("inf")
     maximum = -float("inf")
     for row in matrix:
        try:
            value = float(row[i])
            minimum = min(minimum,value)
            maximum = max(maximum,value)
        except:
            pass
     return minimum,maximum

def removeInvariantColumns(table,header):
    '''
    Remove all columns containing the same value in every row
    of the given table.
    '''  
    invariantColumns        = collections.OrderedDict()
    invariantColumnsIndices = []

    for col in range(0,len(header)):
        current = column(table,col) 
        if all(item.strip()==current[0].strip() for item in current):
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

def createFilterKeysToColumnIndexMapping(filterSet,columnNames):
    """
    Returns the parameter to index mappings for the keys of the given 
    parameterSpace dictionary.
    Performs an one-sided on-the-fly check if all keys are column names.
    """
    indexMapping = collections.OrderedDict()
    
    success = True
    columnNamesSet = set(columnNames)
    for key in filterSet:
        if key in columnNamesSet:
            indexMapping[key] = columnNames.index(key)
        else:
            print("ERROR: parameter key '"+key+"' is not a column name of the table!",file=sys.stderr)
            success = False
    if not success:
        print("ERROR: program aborted since not all parameter keys are a column name of the table.\n",file=sys.stderr)
        print("INFO: found table column names:\n"+"\n".join(columnNames),file=sys.stderr)
        sys.exit()
    
    return indexMapping

def tableFilter(row):
    match = True
    for key,index in filterColumnsToIndices.items():
        match = match and row[index]==rowFilter[key]
    return match

def minMaxFilter(row):
    return float(row[minMaxIndex])==minMax

def getColumnsSortingKey(row):
    keyTuple = ()
    for index in keyIndices:
      try:
          keyTuple += (float(row[index]),)
      except ValueError:
          keyTuple += (row[index],)
    return keyTuple

def parseArgs():
    parser = argparse.ArgumentParser(
        description="This is tableslicer.py: A small tool for extracting columns from a table which is scanned according to a filter.",
            epilog="End of help message.",
    )
    # some mandatory arguments from Toolkit1
    parser.add_argument("--filter", nargs="+", default=["*"],
        help="Specify a list of key-value pairs. Example: ./tableslicer.py --filter order=3 maximumMeshDepth=3 ...")
    parser.add_argument("--cols",   nargs="+",  default=["*"],
        help="Specifiy the list columns you want to read from the rows matching the filter. Example: ./tableslicer.py ... --cols cores realtime_min ")
    
    parser.add_argument("--min", nargs="?", default="",
        help="Specify the column you want to determine the minimum value of. All (filtered) rows with that value will be written out. If you do not specify anything, the last column will be used. Example: ./tableslicer.py --filter .. --cols order maximumMeshDepth ... --min order ")
    
    parser.add_argument("--max", nargs="?", default="",
        help="Specify the column you want to determine the maximum value of. All (filtered) rows with that value will be written out. If you do not specify anything, the last column will be used. Example: ./tableslicer.py --filter .. --cols order maximumMeshDepth ... --max order ")

    parser.add_argument("--header", dest="header", action="store_true",help="Write a header to the output file.")
    parser.add_argument("--no-header", dest="header", action="store_false",help="Write no header to the output file.")
    parser.set_defaults(header=True)
    
    parser.add_argument("--compress", dest="compress", action="store_true",help="Remove columns where the same value is found in every row.")
    parser.add_argument("--no-compress", dest="compress", action="store_false",help="Do not remove columns where the same value is found in every row.")
    parser.set_defaults(compress=False)
    
    parser.add_argument("-s", "--sort", nargs="+", default=[],
        help="Specify a list of sorting key columns. Order is important. Example: ./tableslicer.py ... ... --cols fused cores --sort cores fused ")
    
    parser.add_argument("--input-delim", dest="inputDelim", nargs="?", default=",",
        help="Specify the delimiter used in the input table.")
    
    parser.add_argument("--output-delim", dest="outputDelim", nargs="?", default=",",
        help="Specify the delimiter for the output table.")
    
    parser.add_argument("table",
        type=argparse.FileType("r"),nargs="?",
        help="The CSV table to work with.",
        default=sys.stdin)
    
    parser.add_argument("--output",
        type=argparse.FileType("w"),
        help="The output file.",
        default=sys.stdout)
    
    return parser.parse_args()


if __name__ == "__main__":
    import sys,os
    import csv
    import collections
    import itertools
    import argparse
    
    import math
        
    args = parseArgs()
  
    # read table
    columnNames = next(args.table)
    columnNames = columnNames.strip()
    columnNames = columnNames.split(",")
    tableData   = list(csv.reader(args.table,delimiter=args.inputDelim))
    args.table.close()
   
    ##
    # filter the rows
    ##
    filteredRows = []
    if args.filter[0].strip()!="*":
        # construct row filter
        errorOccured = False
        rowFilter   = {} 
        for item in args.filter:
            tokens = item.split("=")
            if len(tokens)!=2:
                print("ERROR: filter item '{}' does not have shape 'key=value' or 'key=\"value\"'".format(item), file=sys.stderr)
                errorOccured = True
            else:
                rowFilter[tokens[0].strip()]=tokens[1].strip("\"")
        if errorOccured:
            print("Program aborted.", file=sys.stderr)
            sys.exit()
        
        filterColumnsToIndices = createFilterKeysToColumnIndexMapping(rowFilter,columnNames)
        filteredRows = list(filter(tableFilter,tableData))
    else:
        filteredRows = tableData
    
    ##
    # filter the columns
    ##
    if args.cols[0].strip()!="*":
       extractedColumnsToIndices = createFilterKeysToColumnIndexMapping(args.cols,columnNames)
    else:
       extractedColumnsToIndices = createFilterKeysToColumnIndexMapping(columnNames,columnNames)
    header = list(extractedColumnsToIndices.keys())

    result = []
    for row in filteredRows:
        resultRow = []
        for index in extractedColumnsToIndices.values(): # is ordered
            resultRow.append(row[index])
        result.append(resultRow)
  
    ##
    # sort
    ##
    if len(args.sort):
        keyIndices = []
        for item in args.sort:
            if item not in header:
                print("ERROR: sorting key '{}' must be a column of the table: Available columns: {}'".format(item), ", ".join(header))
                sys.exit()
            else:
                keyIndices.append(header.index(item))

        result = sorted(result,key=getColumnsSortingKey)

    ##
    # min,max filter
    ##
    if args.min!="":
        minMaxIndex = -1
        if args.min != None:
           if args.min in extractedColumnsToIndices.keys():
               minMaxIndex = list(extractedColumnsToIndices.keys()).index(args.min)
           else: 
               print("ERROR: column to find minimum in, '{}', must be a column of the sliced table: Available columns: {}'".format(item), ", ".join(extractedColumnsToIndices))
               sys.exit()
        minMax, maximum = minAndMaxInColumn(result,minMaxIndex)
        result = list(filter(minMaxFilter,result))
 
    if args.max!="":
        minMaxIndex = -1
        if args.max != None:
           if args.max in extractedColumnsToIndices.keys():
               minMaxIndex = list(extractedColumnsToIndices.keys()).index(args.max)
           else: 
               print("ERROR: column to find maximum in, '{}', must be a column of the sliced table: Available columns: {}'".format(item), ", ".join(extractedColumnsToIndices))
               sys.exit()
        minimum, minMax = minAndMaxInColumn(result,minMaxIndex)
        result = list(filter(minMaxFilter,result))
         

    ##
    # compress
    ##
    if args.compress:
        result,header,invariantColumns = removeInvariantColumns(result,header)

    csvwriter = csv.writer(args.output,delimiter=args.outputDelim.replace("\\t","\t"))
    if args.header:
        csvwriter.writerow(header)
    csvwriter.writerows(result)
