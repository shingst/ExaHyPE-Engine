"""
.. module:: sweep_analysis
  :platform: Unix, Windows, Mac
  :synopsis: Submodule containing modules for analysing 
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>, 

:synopsis: Generate benchmark suites for ExaHyPE.
"""
import sys
import configparser
import csv
import collections
import re

def parseList(string):
    """
    Decomposes strings like '"val1,val2",val3,"val4,val5"'
    into a list of strings:
    [ 'val1,val2' ,'val3', 'val4,val5' ]
    """
    string = string.replace("\r","").replace("\n","").replace("\t","")
    for line in csv.reader([string],delimiter=","):
      values = line
      for i,value in enumerate(values):
          values[i] = value.replace("\"","")
      return values

def parseEnvironment(config):
    """
    Parse the environment section.
    """
    environmentSpace = {}
    if "environment" in config and len(config["environment"].keys()):
        for key, value in config["environment"].items():
            environmentSpace[key] = parseList(value)
        if "SHAREDMEM" not in environmentSpace:
            print("ERROR: 'SHAREDMEM' missing in section 'environment'.",file=sys.stderr)
            sys.exit()
    else:
        print("ERROR: Section 'environment' is missing or empty! (Must contain at least 'SHAREDMEM'.)",file=sys.stderr)
        sys.exit()
    
    return environmentSpace

def parseParameters(config):
    """
    Parse the parameters section.
    """
    multipleListings = False
    ungroupedParameterSpace = {}
    groupedParameterSpace  = {}
    if config.has_section("parameters"):
        for key, value in config["parameters"].items():
            if key in ungroupedParameterSpace:
               print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
               multipleListings = True
            ungroupedParameterSpace[key] = parseList(value)
        
    if config.has_section("parameters_grouped"):
        for key, value in config["parameters_grouped"].items():
            if key in groupedParameterSpace:
               print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
               multipleListings = True
            if key in ungroupedParameterSpace:
               print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
               multipleListings = True
            groupedParameterSpace[key] = parseList(value) 
    if multipleListings:
        print("ERROR: Some parameters have been found multiple times. Program aborted.",file=sys.stderr)
        sys.exit()

    if len(groupedParameterSpace)==0 and len(ungroupedParameterSpace)==0:
        print("ERROR: Section 'parameters' is missing or empty! (Must contain at least 'dimension' and 'order'.)",file=sys.stderr)
        sys.exit()
    
    parameterSpace = {}
    parameterSpace.update(ungroupedParameterSpace)
    parameterSpace.update(groupedParameterSpace)
    
    # always put None in order to have at least one element
    groupedParameterSpace[None] = [None] 

    blackList = ["tasks","cores","backgroundTasks"]
    for key in blackList:
        if key in parameterSpace:
            print("ERROR: The following keys are reserved: "+",".join(blackList)+".",file=sys.stderr)
            sys.exit()
 
    return ungroupedParameterSpace,groupedParameterSpace,parameterSpace

def parseCoreCounts(coreCounts):
    myCoreCounts = [] 
    for i,cores in enumerate(coreCounts):
        m1 = re.search("([0-9]+)\s*\:\s*([0-9]+)",cores.strip())
        m2 = re.search("([0-9]+)",cores.strip())
        if m1:
            myCoreCounts.append(m1.group(1)+":"+m1.group(2))
        elif m2:
            myCoreCounts.append(m2.group(1)+":"+m2.group(1)) # default
        else:
            print("ERROR: please use either the notation '<cores>' or '<cores>:<background threads>' on the right-hand side of"+\
                 "'cores' or 'cores_grouped'!",file=sys.stderr)
            sys.exit()
    return myCoreCounts


def parseOptionsFile(optionsFile,ignoreMetadata=False):
    configParser = configparser.ConfigParser()
    configParser.optionxform=str
    configParser.read(optionsFile)
    
    general          = dict(configParser["general"])
    exahypeRoot      = general["exahype_root"]
    outputPath       = general["output_path"]
    projectPath      = general["project_path"]
    projectName      = general["project_name"]
    
    buildFolder       = "build"
    buildFolderPath   = exahypeRoot+"/"+outputPath+"/" + buildFolder
    scriptsFolderPath = exahypeRoot+"/"+outputPath+"/scripts"
    resultsFolderPath = exahypeRoot+"/"+outputPath+"/results"
    historyFolderPath = exahypeRoot+"/"+outputPath+"/history"
    
    jobs                                                           = dict(configParser["jobs"])
    environmentSpace                                               = parseEnvironment(configParser)
    ungroupedParameterSpace, groupedParameterSpace, parameterSpace = parseParameters(configParser)

    jobClass   = "unknown"
    islands    = "unknown" 
    if configParser.has_option("jobs","class"):
        jobClass = jobs["class"].strip();
    if configParser.has_option("jobs","islands"):
        islands = jobs["islands"].strip();

    rankCounts = [x.strip() for x in jobs["ranks"].split(",")]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    # cores 
    coreCounts       =["+"]
    coreCountsGrouped=[None]
    if configParser.has_option("jobs","cores"):
        coreCounts = [x.strip() for x in jobs["cores"].split(",")]
        coreCounts = parseCoreCounts(coreCounts)
    elif configParser.has_option("jobs","cores_grouped"):
        coreCountsGrouped = [x.strip() for x in jobs["cores_grouped"].split(",")]
        coreCountsGrouped = parseCoreCounts(coreCountsGrouped)
    if configParser.has_option("jobs","cores") == configParser.has_option("jobs","cores_grouped"):
        print("ERROR: please either specify 'cores' or 'cores_grouped'!",file=sys.stderr)
        sys.exit()
    # runs
    runNumbers       =["+"]
    runNumbersGrouped=[None]
    if configParser.has_option("jobs","run"):
        runNumbers = [x.strip() for x in jobs["run"].split(",")]
    elif configParser.has_option("jobs","run_grouped"):
        runNumbersGrouped = [x.strip() for x in jobs["run_grouped"].split(",")]
    if configParser.has_option("jobs","run") == configParser.has_option("jobs","run_grouped"):
        print("ERROR: please either specify 'run' or 'run_grouped'!",file=sys.stderr)
        sys.exit()
   
    Options = collections.namedtuple("options", \
           ("general jobs "
            "environmentSpace "
            "ungroupedParameterSpace groupedParameterSpace parameterSpace "
            "exahypeRoot outputPath projectPath projectName "
            "buildFolder "
            "buildFolderPath scriptsFolderPath "
            "resultsFolderPath historyFolderPath "
            "jobClass islands "
            "rankCounts nodeCounts "
            "coreCounts coreCountsGrouped "
            "runNumbers runNumbersGrouped"))
    
    options = Options(
      general                 = general,
      jobs                    = jobs,\
      environmentSpace        = environmentSpace,\
      ungroupedParameterSpace = ungroupedParameterSpace,\
      groupedParameterSpace   = groupedParameterSpace,\
      parameterSpace          = parameterSpace,\
      \
      exahypeRoot      = exahypeRoot,\
      outputPath       = outputPath,\
      projectPath      = projectPath,\
      projectName      = projectName,\
      \
      buildFolder       = buildFolder,\
      buildFolderPath   = buildFolderPath,\
      scriptsFolderPath = scriptsFolderPath,\
      resultsFolderPath = resultsFolderPath,\
      historyFolderPath = historyFolderPath,\
      \
      jobClass          = jobClass,\
      islands           = islands,\
      rankCounts        = rankCounts,\
      nodeCounts        = nodeCounts,\
      coreCounts        = coreCounts,\
      coreCountsGrouped = coreCountsGrouped,\
      runNumbers        = runNumbers,\
      runNumbersGrouped = runNumbersGrouped\
    )
    
    return options
