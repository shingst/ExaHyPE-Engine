#!/usr/bin/env python3
def parseArgument(argv,i):
    if i<len(argv):
        return argv[i]
    else:
        return None

def haveToPrintHelpMessage(argv):
    """
    Check if we have to print a help message.
    """
    result = parseArgument(argv,1)==None
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def parseList(string):
    """
    Decomposes strings like '"val1,val2",val3,"val4,val5"'
    into a list of strings:
    [ 'val1,val2' ,'val3', 'val4,val5' ]
    """
    for line in csv.reader([string],delimiter=","):
      values = line
      for i,value in enumerate(values):
          values[i] = value.replace(" ","")
      return values

def parseParameterSpace(config,section):
    """
    Parse the environment section.
    """
    parameterSpace = collections.OrderedDict()
    if section in config and len(config[section].keys()):
        for key, value in config[section].items():
            parameterSpace[key] = parseList(value)
    else:
        print("ERROR: Section '"+section+"' must not be empty!",file=sys.stderr)
        sys.exit()
    
    return parameterSpace

def dictProduct(dictionary):
    """
    Computes the Cartesian product of a dictionary of lists as 
    a list of dictionaries.
    
    Gladly copied this code from:
    https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists
    
    Example input:
    options = {"number": [1,2,3], "color": ["orange","blue"] }
    
    Example output:
    [ {"number": 1, "color": "orange"},
      {"number": 1, "color": "blue"},
      {"number": 2, "color": "orange"},
      {"number": 2, "color": "blue"},
      {"number": 3, "color": "orange"},
      {"number": 3, "color": "blue"}
    ]
    """
    return (collections.OrderedDict(zip(dictionary, x)) for x in itertools.product(*dictionary.values()))

def column(matrix, i):
    return [row[i] for row in matrix]

def autolabel(axes,rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        label = rect.get_label()
        axes.text(rect.get_y() + rect.get_height()/2., 1.05*rect.get_width(),
                "%s" % label, 
                ha='center', va='bottom')

def plot():
    """
    Create a plot per plotDict. 
    Per plot, plot all rows found for the elements of the perPlotSpace.
    """
    for plotDict in dictProduct(plotsSpace):
        # create new plot
        figure = pyplot.figure()
        axes   = figure.add_subplot(111)
        
        plotDictKey = ",".join("%s=%s" %  pair for pair in plotDict.items())
        
        counter = 0
        rects = []
        for perPlotDict in dictProduct(perPlotSpace):
            def tableFilter(row):
                hit = True
                for key,index in indexMappingPlots.items():
                    hit = hit and row[index]==plotDict[key]
                for key,index in indexMappingPerPlot.items():
                    hit = hit and row[index]==perPlotDict[key]
                return hit
            
            filtered = list(filter(tableFilter,tableData))
            perPlotDictKey = ",".join("%s=%s" %  pair for pair in perPlotDict.items())
            filterKey = plotDictKey + perPlotDictKey
            if len(filtered)>1:
                print("ERROR: Parameter combinations are not unique!",file=sys.stderr)
                print("ERROR: Found multiple rows for filter key=("+filterKey+").",file=sys.stderr)
                print("ERROR: Differing column values:",file=sys.stderr)
                for index,name in enumerate(columnNames):
                    values = set(column(filtered,index))
                    if len(values)>1:
                        print("ERROR: "+name+"={"+",".join(values)+"}",file=sys.stderr)
                sys.exit()
            elif len(filtered)==1:
                dataPoint = float(filtered[0][dataColumnIndex])
                rect = axes.barh(-counter,dataPoint,height=0.4,color='0.4',align='center',log=False,label=perPlotDictKey)
                rects.append(rect)
                counter += 1
            elif len(filtered)==0:
                print("WARNING: Found no rows for key=("+filterKey+")!",file=sys.stderr)
        
        if counter>0:
            axes.set_title(plotDictKey.replace(",","   "))
            axes.get_yaxis().set_visible(False)
            axes.set_xlabel(dataColumnName)
            axes.grid(True, which='x')
            
            autolabel(axes,rects)
            
            figure.show()
            input("Press Enter to continue...")

def getDataColumnIndex():
    """
    Returns the parameter to index mappings for the keys of the given 
    parameterSpace dictionary.
    Performs an one-sided on-the-fly check if all keys are column names.
    """
    columnIndex = 0;
    if dataColumnName in columnNames:
        columnIndex = columnNames.index(dataColumnName)
    else:
      print("ERROR: program aborted since data column to plot "+dataColumnName+" is not a column name of the table.",file=sys.stderr)
      print("ERROR: found table column names: "+",".join(columnNames),file=sys.stderr)
      sys.exit()
    
    return columnIndex

def createParameterKeysToColumnIndexMapping(parameterSpace):
    """
    Returns the parameter to index mappings for the keys of the given 
    parameterSpace dictionary.
    Performs an one-sided on-the-fly check if all keys are column names.
    """
    indexMapping   = {}
    
    success = True
    columnNamesSet = set(columnNames)
    for key in parameterSpace:
        if key in columnNamesSet:
            indexMapping[key] = columnNames.index(key)
        else:
            print("ERROR: parameter key '"+key+"' is not a column name of the table!",file=sys.stderr)
            success = False
    if not success:
        print("ERROR: program aborted since not all parameter keys are a column name of the table.",file=sys.stderr)
        print("ERROR: found table column names: "+",".join(columnNames),file=sys.stderr)
        sys.exit()
    
    return indexMapping

if __name__ == "__main__":
    import sys
    import csv
    import configparser
    import collections
    import itertools
    
    import matplotlib
    import matplotlib.pyplot as pyplot
    
    if haveToPrintHelpMessage(sys.argv):
        info = \
"""tableplotter.py:

run:

./tableplotter.py myoptions.ini mytable.csv

NOTE: The order of the parameters in the section 'per_plot' is preserved. 
      You thus have some control over the position of bars in the resulting diagrams.

NOTE: tableplotter preserves the order of the parameters. You thus have control over
      the order of the position of the parameters in the loop.
"""
        print(info) # correctly indented
        sys.exit()
    
    # read options
    optionsFilePath = parseArgument(sys.argv,1)
    tablePath       = parseArgument(sys.argv,2)
    
    configParser = configparser.ConfigParser()
    configParser.optionxform=str
    configParser.read(optionsFilePath)
    
    plotPrefix      = configParser["output"]["prefix"].replace("\"","")
    outputPath      = configParser["output"]["path"].replace("\"","")
    buildFolderPath = outputPath+"/plots"
    
    dataColumnName = configParser["to_plot"]["data"].replace("\"","")
    
    plotsSpace     = parseParameterSpace(configParser,"plots")
    perPlotSpace   = parseParameterSpace(configParser,"per_plot")
    
    # read table
    tableFile   = open(tablePath, 'r')
    columnNames = next(tableFile)
    columnNames = columnNames.strip()
    columnNames = columnNames.split(";")
    tableData   = list(csv.reader(tableFile,delimiter=";"))
    tableFile.close()
    
    dataColumnIndex     = getDataColumnIndex()
    indexMappingPlots   = createParameterKeysToColumnIndexMapping(plotsSpace)
    indexMappingPerPlot = createParameterKeysToColumnIndexMapping(perPlotSpace)
    
    # plot
    plot()
