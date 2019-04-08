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

def printJobTemplate(machine):
    if machine=="supermuc":
        template="""#!/bin/bash
# Mandatory parameters are:
# time, ranks, nodes,
# job_name, output_file, error_file, 
# body
# 
# Optional parameters are:
# tasks, coresPerRank, mail

#@ job_type     = parallel
#@ class        = {{class}}
#@ total_tasks  = {{ranks}}
#@ node         = {{nodes}}
#@ island_count = {{islands}}
#@ network.MPI = sn_all,not_shared,us 
#@ energy_policy_tag = ExaHyPE_Euler_energy_tag
#@ minimize_time_to_solution = yes
#@ wall_clock_limit = {{time}}
#@ job_name = {{job_name}}
#@ error  =  {{error_file}}
#@ output =  {{output_file}}
#@ notification=complete
#@ notify_user={{mail}}
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

export OMP_NUM_THREADS={{coresPerRank}}
export MP_TASK_AFFINITY=core:{{coresPerRank}}

{{body}}"""
        print(template)
    elif machine=="supermuc-skx":
        template="""#!/bin/bash
# Mandatory parameters are:
# time, ranks, nodes,
# job_name, output_file, error_file, 
# body
# 
# Optional parameters are:
# tasks, coresPerRank, mail

# see: https://doku.lrz.de/display/PUBLIC/Job+Processing+with+SLURM+on+SuperMUC-NG

#SBATCH --job-name={{job_name}}
#SBATCH -o {{output_file}}
#SBATCH -e {{error_file}}
#SBATCH -t {{time}}
#SBATCH --partition={{class}}
#SBATCH --ntasks={{ranks}}
#SBATCH --nodes={{nodes}}
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --mail-user={{mail}}
#SBATCH --mail-type=END
#SBATCH --export=NONE
#SBATCH --account=pr48ma
#SBATCH --chdir ./

# due to the fact that SLURM performs settings of certain variables in its prologue that can cause MPI crashes, 
# and cannot be overridden in the default module mechanism, the following module should be loaded in SLURM scripts before executing parallel programs:
module load slurm_setup

module load intel/19.0
module load gcc/7
module load tbb/2019
module load mpi.intel/2019

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

{{body}}"""
        print(template)
    elif machine=="hamilton":
        template="""#!/bin/bash
# Mandatory parameters are:
# time, ranks, nodes,
# job_name, output_file, error_file, 
# body
# 
# Optional parameters are:
# tasks, coresPerRank, mail

#SBATCH --job-name={{job_name}}
#SBATCH -o {{output_file}}
#SBATCH -e {{error_file}}
#SBATCH -t {{time}}
#SBATCH --exclusive
#SBATCH -p par7.q
#SBATCH --mem=MaxMemPerNode
#SBATCH --ntasks={{ranks}}
#SBATCH --nodes={{nodes}}
#SBATCH --mail-user={{mail}}
#SBATCH --mail-type=END
module purge
module load slurm
module load intelmpi/intel/2018.2
module load intel/xe_2018.2
module load gcc/7.3.0
module load python/3.6.8
module load likwid

export TBB_SHLIB="-L/ddn/apps/Cluster-Apps/intel/xe_2018.2/tbb/lib/intel64/gcc4.7 -ltbb"

export I_MPI_FABRICS="tmi"

{{body}}"""
        print(template)
    else:
        print("ERROR: No job template found for machine {}. Available options: '{}'".format(machine,"','".join(["supermuc","supermuc-skx","hamilton"])),file=sys.stderr)

def printOptionsFileTemplate():
    template="""[general]
exahype_root   = EXAHYPE_ROOT
project_name   = PROJECT
project_path   = PROJECT_PATH

spec_template    = %(project_path)s/SPEC_TEMPLATE.exahype2-template 
job_template     = %(project_path)s/JOB_TEMPLATE.job-template

output_path      = %(project_path)s/OUTPUT_PATH
make_threads     = MAKE_THREADS

run_command      = mpiexec|poe

likwid           = LIKWID<OPTIONAL>

job_submission   = llsubmit|sbatch
job_cancellation = llcancel|scancel

compile_time_parameters = COMPILE_TIME_PARAMETERS

[jobs]
time = TIME
mail = EMAIL

num_cpus = NUM_CPUS

class             = CLASS<OPTIONAL>
islands           = ISLANDS<OPTIONAL>

ranks_nodes_cores = RANKSxNODESx{CORES:BACKGROUND_THREADS},
                    1x1x{1:1},
                    ...

run  = RUN<OPTION1>
run_grouped = RUN_GROUPED<OPTION2>

[environment]
;<ALL OPTIONAL>
EXAHYPE_CC      = mpiCC
EXAHYPE_FC      = mpif90
COMPILER        = Intel
MODE            = RELEASE
DISTRIBUTEDMEM  = MPI
SHAREDMEM       = TBB
USE_IPO         = on
COMPILER_CFLAGS = ""

[parameters]
;<ALL BELOW MAY BE MOVED TO parameters_grouped or vice versaL>
architecture = ARCHITECTURE
dimension    = DIMENSION
;<ALL BELOW OPTIONAL>
param0       = PARAM0

[parameters_grouped]
;<ALL BELOW OPTIONAL>
param0 = PARAM0"""
    print(template)
