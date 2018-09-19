#!/usr/bin/env python3

def parseArgs():
    parser = argparse.ArgumentParser(
        description="This is speedupcalculatorpy: A small tool for computing speedups for data stored in a CSV table. The last column of a table row is assumed to be the data column. The first row is assumed to store reference value.",
            epilog="End of help message.",
    )
    parser.add_argument('--header', dest='header', action='store_true',help="Write a header to the output file.")
    parser.add_argument('--no-header', dest='header', action='store_false',help="Write no header to the output file.")
    parser.set_defaults(header=True)
    
    parser.add_argument("--enumerate-keys", dest="enumerateKeys", action="store_true",help="Enumerate the key columns instead of including them.")
    parser.add_argument("--no-enumerate-keys", dest="enumerateKeys", action="store_false",help="Include the key columns.")
    parser.set_defaults(enumerateKeys=False)
    
    parser.add_argument("--data", dest="keepData", action="store_true",help="Include the original data column into the output file.")
    parser.add_argument("--no-data", dest="keepData", action="store_false")
    parser.set_defaults(keepData=True)
    
    parser.add_argument("-t", "--tikz", dest="forTikz", action="store_true",help="Has the effect of --enumerate-keys --no-headers --no-data'. Overwrites these options.")
    parser.set_defaults(forTikz=False)
    
    parser.add_argument("table",
       type=argparse.FileType("r"),nargs="?",
       help="The CSV table to work on.",
       default=sys.stdin)
    
    parser.add_argument("-o","--output",
       type=argparse.FileType("w"),nargs="?",
       help="Output file",
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

    enumerateKeys = args.enumerateKeys
    header        = args.header
    keepData      = args.keepData
    forTikz       = args.forTikz
    table         = args.table
    output        = args.output

    if forTikz:
       enumerateKeys = True
       header        = False
       keepData      = False
 
    # read table
    firstRow  = next(args.table)
    firstRow  = firstRow.strip().split(",")
    tableData = list(csv.reader(args.table,delimiter=","))
    args.table.close()

    # check if has header; include header into output if requested 
    result = []
    try:
        float(firstRow[-1])
        tableData.insert(0,firstRow)
    except:
       if header:
           if not enumerateKeys:
               if not keepData:
                  result.append(firstRow[:-1])
               else:
                  result.append(firstRow)
               result[0].append("speedup")
           else:
               if not keepData:
                   result.append(["n","speedup"])
               else:
                   result.append(["n",firstRow[-1],"speedup"])
 
    
    reference = float(tableData[0][-1])
    n = 0
    for row in tableData:
        resultRow = None
        if not enumerateKeys:
            resultRow = row
        else:
            resultRow = [n,row[-1]]
        if not keepData:
            resultRow = resultRow[:-1]
        resultRow.append("%1.2f" % ( reference/float(row[-1])) ) # compute speedup
        result.append(resultRow)
        n += 1
    
    csvwriter = csv.writer(output)
    csvwriter.writerows(result)
