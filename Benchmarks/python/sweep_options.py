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
    if "parameters" in config and len(config["parameters"].keys()):
        for key, value in config["parameters"].items():
            ungroupedParameterSpace[key] = parseList(value)
            if key in ungroupedParameterSpace:
               print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
               multipleListings = True
        
        if config.has_section("parameters_grouped"):
            for key, value in config["parameters_grouped"].items():
                groupedParameterSpace[key] = parseList(value)
                if key in groupedParameterSpace:
                   print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
                   multipleListings = True
                if key in ungroupedParameterSpace:
                   print("ERROR: Parameter '"+key+"' found multiple times.",file=sys.stderr)
                   multipleListings = True
        
        if multipleListings:
            print("ERROR: Some parameters have been found multiple times. Program aborted.",file=sys.stderr)
            sys.exit()
        if "order" not in parameterSpace:
            print("ERROR: 'order' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        elif "dimension" not in parameterSpace:
            print("ERROR: 'dimension' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        elif "optimisation" not in parameterSpace:
            print("ERROR: 'optimisation' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        elif "architecture" not in parameterSpace:
            print("ERROR: 'architecture' missing in section 'parameters'.",file=sys.stderr)
            sys.exit()
        
        # if there are no grouped parameters specified put a None key into the dictionary
        if len(groupedParameterSpace)==0:
            groupedParameterSpace[None] = []
            
        parameterSpace = {}
        parameterSpace.update(groupedParameterSpace)
        parameterSpace.update(groupedParameterSpace)
    else:
        print("ERROR: Section 'parameters' is missing or empty! (Must contain at least 'dimension' and 'order'.)",file=sys.stderr)
        sys.exit()
    
    return ungroupedParameterSpace,groupedParameters,parameterSpace

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
    
    jobs             = dict(configParser["jobs"])
    environmentSpace                                               = parseEnvironment(configParser)
    ungroupedParameterSpace, groupedParameterSpace, parameterSpace = parseParameters(configParser)

    jobClass   = "unknown"
    islands    = "unknown" 
    if configParser.has_option("jobs","class"):
        jobClass = jobs["class"].strip();
    if configParser.has_option("jobs","islands"):
        islands = jobs["islands"].strip();
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = [x.strip() for x in jobs["cores"].split(",")]
    runNumbers = [x.strip() for x in jobs["run"].split(",")]
    
    Options = collections.namedtuple("options", \
           ("general jobs environmentSpace "
            "ungroupedParameterSpace groupedParameterSpace parameterSpace "
            "exahypeRoot outputPath projectPath projectName "
            "buildFolder "
            "buildFolderPath scriptsFolderPath "
            "resultsFolderPath historyFolderPath "
            "jobClass islands "
            "nodeCounts taskCounts coreCounts runNumbers"))
    
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
      jobClass   = jobClass,\
      islands    = islands,\
      nodeCounts = nodeCounts,\
      taskCounts = taskCounts,\
      coreCounts = coreCounts,\
      runNumbers = runNumbers\
    )
    
    return options
