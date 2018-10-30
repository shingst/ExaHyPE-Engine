"""
.. module:: sweep analysis
  : platform: Unix, Windows, Mac
  : synopsis: Submodule containing modules for analysing

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
          values[i] = value.strip().replace("\"","")
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

    blackList = ["ranksPerNode","coresPerRank","consumerTasks"]
    for key in blackList:
        if key in parameterSpace:
            print("ERROR: The following keys are reserved: "+",".join(blackList)+".",file=sys.stderr)
            sys.exit()
    
    # compile-time parameters
    compileTimeParameterSpace = collections.OrderedDict()
    print (config.has_option("general","compile_time_parameters"))
    if  config.has_option("general","compile_time_parameters"):
        compileTimeParameters = parseList(config.get("general", "compile_time_parameters"))
        for key in compileTimeParameters:
            if key in parameterSpace:
                if key not in ["dimension","architecture"]:
                    compileTimeParameterSpace[key] = parameterSpace[key]
            else:
                print("ERROR: Did not find 'compile_time_parameters' entry '{}' in section 'parameters' or 'grouped_parameters'.".format(key),file=sys.stderr)
                sys.exit()
    else:
        print("ERROR: Did not find required 'compile_time_parameters' option where you specify the parameters whose variation requires a rerun of the toolkit and a recompilaton of the project folder sources!" +\
              "Note that you do not need to specify 'dimensions' and 'architecture' as they require a rebuild of Peano as well. They are treated separately.",file=sys.stderr)
        sys.exit()
 
    return ungroupedParameterSpace,groupedParameterSpace,parameterSpace,compileTimeParameterSpace

def compareRanksNodesCoreCountsWithEachOther(ranksNodesCoreCountsList):
    overlap = False
    for i,tuple1 in enumerate(ranksNodesCoreCountsList):
        for j,tuple2 in enumerate(ranksNodesCoreCountsList):
            if i != j:
                coresNodesIdentical = tuple1.nodes==tuple2.nodes and\
                                      tuple1.ranks==tuple2.ranks
                if coresNodesIdentical:
                    for coreCount1 in tuple1.coreCounts:
                        for coreCount2 in tuple2.coreCounts:
                            if coreCount1.cores==coreCount2.cores and\
                               coreCount1.consumers==coreCount2.consumers:
                                 print("ERROR: Detected {<cores>:<consumers>} overlap for the two 'ranks_nodes_cores' entries '"+tuple1.text+"' and '"+tuple2.text+"'.",file=sys.stderr)
                                 overlap = True;
    return overlap

def parseRanksNodesCoreCounts(jobs):
    """
    Returns a list of tuples of the form (ranks,nodes,coreCounts).
    coreCounts is an ordered dict with key 'cores' and value 'background job'
    consumer tasks.
    """
    RanksNodesCoreCounts = collections.namedtuple("ranksNodesCoreCounts", \
        "nodes ranks coreCounts text")
    CoreCount = collections.namedtuple("coreCount", \
        "cores consumers")
    # example: "29 x 29 x {4:2,8:4},\\n758 x 29 x {1:1}"
    entryPattern      = re.compile(r"(\s*([0-9]+)\s*x\s*([0-9]+)\s*x\s*\{(([0-9]|\:|,|\s)+)\})");
    # example: 4:2,8:4
    # example: 4:2
    coreCountsPattern = re.compile(r"([0-9]+)\s*\:\s*([0-9]+)");
    
    ranksNodesCoreCountsList = []
    for entryMatch in re.findall(entryPattern, jobs["ranks_nodes_cores"]):
        coreCountsList = []
        for coresMatch in re.findall(coreCountsPattern,entryMatch[3]):
            coreCountsList.append( CoreCount( cores=coresMatch[0], consumers=coresMatch[1] ) )
        if not coreCountsList:
            print("ERROR: the cores specification (the part in curly braces {...}) in 'ranks_nodes_cores' does not have the form: '{<cores0>:<consumertasks0>,<cores1>:<consumertasks1>',...}'.\n"\
                  "       Valid examples: '{1:1}', '{1:1,2:2}'",file=sys.stderr)
            sys.exit()
        ranks = entryMatch[1]
        nodes = entryMatch[2]
        if int(nodes) > int(ranks):
            print("ERROR: specified ranks ("+ranks+") must always be greater than or equal to specified nodes ("+nodes+").\n"\
                  "       Error occured in entry: '"+entryMatch[0]+"'",file=sys.stderr) 
            sys.exit()
        ranksNodesCoreCountsList.append( RanksNodesCoreCounts( ranks=ranks, nodes=nodes, coreCounts=coreCountsList, text=entryMatch[0] ) )
        
    if not ranksNodesCoreCountsList:
        print("ERROR: could not find any 'ranks_nodes_cores' entries in the form: 'ranks_nodes_cores = <ranks> x <nodes> x {<cores0>:<consumertasks0>,<cores1>:<consumertasks1>,...}'.\n"\
              "       Valid examples: 'ranks_nodes_cores = 29 x 29 x {1:1}', 'ranks_nodes_cores = 758 x 29 x {8:4,16:8}'",file=sys.stderr)
        sys.exit()
    elif compareRanksNodesCoreCountsWithEachOther(ranksNodesCoreCountsList):
        sys.exit()
    return ranksNodesCoreCountsList

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
    
    jobs = dict(configParser["jobs"])
    environmentSpace = parseEnvironment(configParser)
    ungroupedParameterSpace, groupedParameterSpace, parameterSpace, compileTimeParameterSpace = parseParameters(configParser)

    jobClass   = "unknown"
    islands    = "unknown" 
    if configParser.has_option("jobs","class"):
        jobClass = jobs["class"].strip();
    if configParser.has_option("jobs","islands"):
        islands = jobs["islands"].strip();
    
    ranksNodesCoreCounts = parseRanksNodesCoreCounts(jobs)
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
            "compileTimeParameterSpace "
            "exahypeRoot outputPath projectPath projectName "
            "buildFolder "
            "buildFolderPath scriptsFolderPath "
            "resultsFolderPath historyFolderPath "
            "jobClass islands "
            "ranksNodesCoreCounts "
            "runNumbers runNumbersGrouped"))
    
    options = Options(
      general                   = general,
      jobs                      = jobs,\
      environmentSpace          = environmentSpace,\
      ungroupedParameterSpace   = ungroupedParameterSpace,\
      groupedParameterSpace     = groupedParameterSpace,\
      parameterSpace            = parameterSpace,\
      compileTimeParameterSpace = compileTimeParameterSpace,\
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
      jobClass             = jobClass,\
      islands              = islands,\
      ranksNodesCoreCounts = ranksNodesCoreCounts,\
      runNumbers           = runNumbers,\
      runNumbersGrouped    = runNumbersGrouped\
    )
    
    return options

