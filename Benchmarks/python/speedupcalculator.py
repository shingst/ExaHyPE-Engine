#!/usr/bin/env python3

def parseArgs():
    parser = argparse.ArgumentParser(
        description="This is speedupcalculatorpy: A small tool for computing speedups for data stored in a CSV table. The last column of a table row is assumed to be the data column. The first row is assumed to store reference value.",
            epilog="End of help message.",
    )
    parser.add_argument('--header', dest='header', action='store_true',help="Write a header to the output file.")
    parser.add_argument('--no-header', dest='header', action='store_false',help="Write no header to the output file.")
    parser.set_defaults(header=True)
    
    parser.add_argument("--keys", dest="keepKeys", action="store_true",help="Include the key columns to the output.")
    parser.add_argument("--no-keys", dest="keepKeys", action="store_false",help="Remove the key columns from the output.")
    parser.set_defaults(keepKeys=True)
    
    parser.add_argument("--multiply-keys", dest="multiplyKeys", action="store_true",help="Include the multiplied value of the keys as a column.")
    parser.add_argument("--no-multiply-keys", dest="multiplyKeys", action="store_false",help="Do not include the multiplied value of the keys as a column.")
    parser.set_defaults(multiplyKeys=False)
    
    parser.add_argument("--data", dest="keepData", action="store_true",help="Include the original data column into the output file.")
    parser.add_argument("--no-data", dest="keepData", action="store_false",help="Remove the original data column from the output file.")
    parser.set_defaults(keepData=True)
    
    parser.add_argument("-t", "--tikz", dest="forTikz", action="store_true",help="Has the effect of --keys --no-headers --no-data'. Overwrites these options.")
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

    header       = args.header
    keepKeys     = args.keepKeys
    multiplyKeys = args.multiplyKeys
    keepData     = args.keepData
    table        = args.table
    output       = args.output
    
    forTikz      = args.forTikz

    if forTikz:
       keepKeys     = False
       multiplyKeys = True
       header       = False
       keepData     = False
 
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
           if keepKeys:
               result.append(firstRow[:-1])
           if multiplyKeys:
               result[0].append("productOfKeys")
           if keepData:
               result[0].append(firstRow[-1])
           result[0].append("speedup")
 
   
    if len(tableData) and len(tableData[0]):
        reference = float(tableData[0][-1])
        n = 0
        for row in tableData:
            resultRow = []
            if keepKeys:
                resultRow += row[:-1]
            if multiplyKeys:
                product = 1
                for key in row[:-1]:
                    try:
                        product *= int(key) 
                    except:
                        print("ERROR: Cannot multiply key column entries as '{}' is not an integer".format(key), file=sys.sterr)
                resultRow.append(str(product))
            if keepData:
                resultRow.append(row[-1])
            
            resultRow.append("%1.2f" % ( reference/float(row[-1])) ) # compute speedup
            result.append(resultRow)
            n += 1
    
    csvwriter = csv.writer(output)
    csvwriter.writerows(result)
