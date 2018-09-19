#!/usr/bin/env python3
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
        match = match and row[index].startswith(rowFilter[key])
    return match

def parseArgs():
    parser = argparse.ArgumentParser(
        description="This is tableslicer.py: A small tool for extracting columns from a table which is scanned according to a filter.",
            epilog="End of help message.",
    )
    # some mandatory arguments from Toolkit1
    parser.add_argument("-f", "--filter", nargs="+", default=["*"],
        help="Specify a list of key-value pairs. Example: ./tableslicer.py --filter order=3 maximumMeshDepth=3 ...")
    parser.add_argument("-c", "--cols",   nargs="+",  default=["*"],
        help="Specifiy the list columns you want to read from the rows matching the filter. Example: ./tableslicer.py ... --cols cores realtime_min ")
    
    parser.add_argument('--header', dest='header', action='store_true',help="Write a header to the output file.")
    parser.add_argument('--no-header', dest='header', action='store_false',help="Write no header to the output file.")
    parser.set_defaults(header=True)
    
    parser.add_argument('table',
        type=argparse.FileType('r'),
        help="The CSV table to work on")
    
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
    tableData   = list(csv.reader(args.table,delimiter=","))
    args.table.close()
   
    # filter the rows
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
    
    # filter the columns
    if args.cols[0].strip()!="*":
       extractedColumnsToIndices = createFilterKeysToColumnIndexMapping(args.cols,columnNames)
    else:
       extractedColumnsToIndices = createFilterKeysToColumnIndexMapping(columnNames,columnNames)

    result = []
    if args.header:
       result.append(list(extractedColumnsToIndices.keys()))
    for row in filteredRows:
        resultRow = []
        for index in extractedColumnsToIndices.values(): # is ordered
            resultRow.append(row[index])
        result.append(resultRow)
    
    csvwriter = csv.writer(sys.stdout)
    csvwriter.writerows(result)
