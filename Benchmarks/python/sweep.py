#!/usr/bin/env python3
"""
.. module:: sweep
  :platform: Unix, Windows, Mac
  :synopsis: Generate benchmark suites for ExaHyPE.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>, 

:synopsis: Generate benchmark suites for ExaHyPE.
"""
def parseArgument(argv,i):
    if i<len(argv):
        return argv[i]
    else:
        return None

def haveToPrintHelpMessage(argv):
    """
    Check if we have to print a help message.
    """
    result = parseArgument(argv,2) not in subprograms or \
             parseArgument(argv,1)==None
    result = result and parseArgument(argv,1) not in ["iniTemplate","jobTemplate"]
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def dictProduct(dicts):
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
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

def hashDictionary(dictionary):
    """
    Hash a dictionary. Sort the dictionary according
    to the keys beforehand.
    """
    chain = ""
    for key,value in sorted(dictionary.items()):
        chain += key+","+value+";"
    
    return hashlib.md5(chain.encode()).hexdigest()[0:8]

def clean(subFolder=""):
    """
    Clean the complete output folder or just a subfolder
    if specified.
    """
    folder = exahypeRoot+"/"+outputPath+"/"+subFolder
    print("rm -r "+folder)
    subprocess.call("rm -r "+folder, shell=True)


def renderSpecFile(templateBody,parameterDict,ranksPerNode,coresPerRank,consumerTasks):
    renderedFile = templateBody
    
    context = dict(parameterDict)
    
    consistent = True
    # verify mandatory options file parameters can be found in template file
    global createdFirstSpecFile
    if not createdFirstSpecFile:
        keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",templateBody)]
        for key in context:
            if key not in keysInTemplate:
                consistent = False
                print("ERROR: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr) 
    # optional parameters
    context["ranksPerNode"]  = ranksPerNode
    context["coresPerRank"]  = coresPerRank
    context["consumerTasks"] = consumerTasks

    if not createdFirstSpecFile:
        for key in context:
            if key not in keysInTemplate:
                print("WARNING: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr)
        # verify template parameters are defined in options file
        for key in keysInTemplate:
            if key not in context:
                consistent = False
                print("ERROR: specification file template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
        if not consistent:
            print("ERROR: subprogram aborted as specification file template and sweep options file are inconsistent.",file=sys.stderr)
            sys.exit("Inconsistent Specification and Sweep options")
        createdFirstSpecFile = True
    
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def verifyLogFilterExists(justWarn=False):
    foundLogFilter = False
    for file in os.listdir(exahypeRoot + "/" + projectPath):
        foundLogFilter = foundLogFilter or file.endswith(".log-filter")

    messageType = "ERROR"
    if justWarn:
        messageTypeV = "WARNING"
    if not foundLogFilter:
        print(messageType+": no 'exahype.log-filter' file could be found in the project folder",file=sys.stderr)
        if not justWarn:
            sys.exit()

def verifyEnvironmentIsCorrect(justWarn=False):
    environmentIsCorrect = True
    
    messageType = "ERROR"
    if justWarn:
        messageType = "WARNING"
    
    cpus = jobs["num_cpus"]
    for environmentDict in dictProduct(environmentSpace):
        for key,value in environmentDict.items():
            os.environ[key]=value
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            if (os.environ["DISTRIBUTEDMEM"].strip() not in ["MPI"]) and int(ranks)>1:
                print(messageType+": DISTRIBUTEDMEM environment variable set to "+environmentDict["DISTRIBUTEDMEM"]+" and ranks is set to "+ranks+" > 1",file=sys.stderr)
                environmentIsCorrect = False
            for coreCount in config.coreCounts:
                cores = coreCount.cores
                if (os.environ["SHAREDMEM"].strip() not in ["TBB","CPP14","OMP","TBBInvade"]) and int(cores)>1:
                    print(messageType+": SHAREDMEM environment variable set to "+environmentDict["SHAREDMEM"]+" and cores set to "+cores+" > 1",file=sys.stderr)
                    environmentIsCorrect = False
                        
    if not justWarn and not environmentIsCorrect:
        print("ERROR: subprogram failed as environment variables are not chosen setup correctly. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the subprogram.",file=sys.stderr)
        sys.exit("Environment setup not chosen correctly")

def verifyAllRequiredParametersAreGiven(specFileTemplate):
    if "dimension" not in parameterSpace:
        print("ERROR: 'dimension' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
        sys.exit("Dimensions not found")
    elif "architecture" not in parameterSpace:
        print("ERROR: 'architecture' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
        sys.exit("Architecture not found")

def unlink():
    """
    Remove old symlinks from the build directory.
    """
    if os.path.exists(buildFolderPath):
        print("clean build directory")
        for file in os.listdir(buildFolderPath):
            if os.path.islink(buildFolderPath+"/"+file):
               os.unlink(buildFolderPath+"/"+file)
        print("")

def link():
    """
    Add user specified symlinks to the build directory.
    """ 
    if os.path.exists(buildFolderPath):
        # symlink local files into build folder
        whitelistFileKeyName = "runtime_dependencies"
        if whitelistFileKeyName in general:
            whiteListFiles = sweep_options.parseList(general[whitelistFileKeyName])
            if len(whiteListFiles) < 1:
                print("[WARNING]    runtime_dependencies_file appears to be empty")
            for f in whiteListFiles:
                fStripped  = f.strip()
                source = None
                if fStripped.startswith("/") or fStripped.startswith("~"):
                    source = fStripped
                else:
                    source = os.path.join(exahypeRoot+"/"+projectPath+"/"+fStripped)
                
                if os.path.isfile(source):
                    dest = os.path.join(buildFolderPath,fStripped)
                    os.symlink(source, dest)
                    print("Symlinked {} to {}".format(source, dest))
                else:
                    print("[ERROR]    could not find source file {} to create symlink".format(source))
                    sys.exit(1)

def getExecutableName(environmentDictHash,architecture,dimension,parameterDict):
    """
    Derives an executable name.
    """
    suffix = architecture+"-d" + dimension

    for key in compileTimeParameterSpace:
        suffix += "-{}_{}".format(key,parameterDict[key])
    
    return buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+suffix

def build(buildOnlyMissing=False, skipMakeClean=False):
    """
    Build the executables.
    
    Args:
    buildOnlyMissing(bool):
       Build only missing executables.
    
    skipMakeClean(bool):
        Do not run make clean, only clean application folder
    """
    subprocess.call("module list",shell=True)
    print("")
    print("found preset ExaHyPE environment variables:")
    exahypeEnv = ["COMPILER", "MODE", "SHAREDMEM", "DISTRIBUTEDMEM", "EXAHYPE_CC", "PROJECT_C_FLAGS", "PROJECT_L_FLAGS", "COMPILER_CFLAGS", "COMPILER_LFLAGS", "FCOMPILER_CFLAGS", "FCOMPILER_LFLAGS"]
    for variable in exahypeEnv:
        if variable in os.environ:
            print(variable+"="+os.environ[variable])
    print("")

    if not os.path.exists(buildFolderPath):
        print("create build directory "+buildFolderPath)
        os.makedirs(buildFolderPath)
    print("")

    unlink()
 
    oldExecutable = exahypeRoot + "/" + projectPath+"/ExaHyPE-"+projectName
    
    verifyLogFilterExists(justWarn=True)        
    verifyEnvironmentIsCorrect(justWarn=True)
    
    architectures = parameterSpace["architecture"]
    dimensions    = parameterSpace["dimension"]
    
    # we will replace values in this dict with compile time parameters
    buildParameterDict = list(dictProduct(parameterSpace))[0]
        
    firstIteration = True
    executables=0
    for environmentDict in dictProduct(environmentSpace):
        for key,value in environmentDict.items():
            os.environ[key]=value
        environmentDictHash = hashDictionary(environmentDict)
        
        for architecture in architectures:
            for dimension in dimensions:
                if not firstIteration and not skipMakeClean:
                    command = "make clean"
                    print(command)
                    process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                 
                for compileTimeParameterDict in dictProduct(compileTimeParameterSpace):
                    executable = getExecutableName(environmentDictHash,architecture,dimension,compileTimeParameterDict)
                    
                    if not os.path.exists(executable) or not buildOnlyMissing:
                        buildParameterDict["architecture"] = architecture
                        buildParameterDict["dimension"]    = dimension
                        for key,value in compileTimeParameterDict.items():
                            buildParameterDict[key] = value
                        
                        buildSpecFileBody = renderSpecFile(specFileTemplate,buildParameterDict,"1","1","1")
                        
                        # build spec file name    
                        suffix = architecture+"-d" + dimension
                        for key in compileTimeParameterSpace:
                            suffix += "-{}_{}".format(key,buildParameterDict[key])
                        buildSpecFilePath = outputPath+"/"+buildFolder+"/"+projectName+"-"+suffix+".exahype"
                            
                        with open(exahypeRoot + "/" + buildSpecFilePath, "w") as buildspecFile:
                            buildspecFile.write(buildSpecFileBody)
                            
                        infoMessage = "building executable with: \n- environment={}\n- architecture={}\n- dimension={}".format(str(environmentDict),architecture,dimension)
                        for key,value in compileTimeParameterDict.items():
                            infoMessage += "\n- {}={}".format(key,value)
                        print(infoMessage)

                        
                        # clean application folder only
                        command = "rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels"
                        print(command)
                        process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        process.communicate()
                        process.wait()
                       
                        # run toolkit
                        toolkitCommand = "{0}/Toolkit/toolkit.sh --format=json -s {0}/{1}".format(exahypeRoot,buildSpecFilePath)
                        print(toolkitCommand,end="",flush=True)
                        process = subprocess.Popen([toolkitCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, toolkitErr) = process.communicate()
                        process.wait()
                        print(output)
                        if toolkitErr: # do not use -v as info output is piped into stderr
                            print(" [FAILED]")
                            print("toolkit output=\n"+output.decode('UTF-8'),file=sys.stderr)
                            print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                            sys.exit("Toolkit failed")
                        
                        if firstIteration and not skipMakeClean:
                            command = "make clean"
                            print(command)
                            process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            (output, err) = process.communicate()
                            process.wait()
                        
                        # call make
                        make_threads=general["make_threads"]
                        makeCommand="make -j"+make_threads
                        print(makeCommand,end="",flush=True)
                        process = subprocess.Popen([makeCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, makeErr) = process.communicate()
                        process.wait()
                        if "build of ExaHyPE successful" in str(output):
                            print(" [OK]")
                        else:
                            print(" [FAILED]",file=sys.stderr)
                            print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                            sys.exit("Build failed")

                        if not os.path.exists(oldExecutable):
                            print("ERROR: could not find built executable '"+oldExecutable+"'. The parameter 'project' in your configuration file is probably wrong." ,file=sys.stderr)
                            sys.exit("Build failed")
                        try:
                            os.rename(oldExecutable,executable)
                        except OSError:
                            import shutil
                            shutil.copy(oldExecutable,executable)
                            os.remove(oldExecutable)
                        
                        print("created executable:"+executable)
                         
                        print("SUCCESS!")
                        print("--------------------------------------------------------------------------------")
                        print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                        print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                        print("--------------------------------------------------------------------------------")
                        executables+=1
                        
                        firstIteration = False
                    else:
                        print("skipped building of '"+executable+"' as it already exists.")

    link()   
 
    print("built executables: "+str(executables))


def renderJobScript(jobScriptTemplate,jobScriptBody,jobs,
                    jobName,jobScriptFilePath,outputFileName,errorFileName,
                    ranks,nodes,ranksPerNode,cores): # cores still necessary?
    """
    Render a job script.
    """
    renderedFile = jobScriptTemplate
    
    context = {}
    # mandatory
    context["ranks"]       = ranks
    context["nodes"]       = nodes
    context["output_file"] = outputFileName
    context["error_file"]  = errorFileName
    context["job_name"]    = jobName 
    
    context["body"] = jobScriptBody   
    
    consistent = True
    # verify all mandatory(!) sweep options are defined in template
    global createdFirstJobScript
    if not createdFirstJobScript:
        keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",jobScriptTemplate)]
        for key in context:
            if key not in keysInTemplate:
                consistent = False
                print("ERROR: parameter '{{"+key+"}}' not found in job script template!",file=sys.stderr)
    
    # put optional sweep options in context
    context["mail"]         = jobs["mail"]
    context["ranksPerNode"]        = ranksPerNode
    context["time"]         = jobs["time"]
    context["class"]        = jobClass
    context["islands"]      = islands
    context["coresPerRank"] = str( int ( int(jobs["num_cpus"]) / int(ranksPerNode) ) )
    
    # now verify template parameters are defined in options file
    if not createdFirstJobScript:
        for key in keysInTemplate:
            if key not in context:
                consistent = False
                print("ERROR: job script template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
        if not consistent:
            print("ERROR: subprogram aborted as job script template and sweep options file are inconsistent.",file=sys.stderr)
            sys.exit("Inconsistent Specification and Sweep options")
    
    context["body"] = jobScriptBody 
 
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def verifyAllExecutablesExist(justWarn=False):
    """
    Verify that all executables exist.
    """    
    architectures = parameterSpace["architecture"]
    dimensions    = parameterSpace["dimension"]
    
    messageType = "ERROR"
    if justWarn:
      messageType = "WARNING"
    
    if not justWarn and not os.path.exists(buildFolderPath):
        print("ERROR: build folder '"+buildFolderPath+"' doesn't exist! Please run subprogram 'build' beforehand.",file=sys.stderr)
        sys.exit("Build folder doesn't exists'")
    
    allExecutablesExist = True
    for environmentDict in dictProduct(environmentSpace):
        environmentDictHash = hashDictionary(environmentDict)
        for architecture in architectures:
            for dimension in dimensions:
                for compileTimeParameterDict in dictProduct(compileTimeParameterSpace):
                    executable = getExecutableName(environmentDictHash,architecture,dimension,compileTimeParameterDict)
                    
                    if not os.path.exists(executable):
                        allExecutablesExist = False
                        infoMessage = "application for the following settings does not exist: \n- environment={}\n- architecture={}\n- dimension={}".format(str(environmentDict),architecture,dimension)
                        for key,value in compileTimeParameterDict.items():
                            infoMessage += "\n- {}={}".format(key,value)
                        print("{}: {}\n('{}')".format(messageType,infoMessage,executable),file=sys.stderr)
            
    if not justWarn and not allExecutablesExist:
        print("ERROR: subprogram failed as not all executables exist. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the 'build' subprogram.",file=sys.stderr)
        sys.exit("Not all executables exists")

def verifySweepAgreesWithHistoricalExperiments():
    """
    If there are any previous experiments ensure that the sweep 
    parameter spaces contain the same axes.
    """
    if os.path.exists(historyFolderPath):
        previousSweeps = [f for f in os.listdir(historyFolderPath) if f.endswith(".ini")]
        print(previousSweeps)
    
        for f in previousSweeps:
            otherOptions = sweep_options.parseOptionsFile(historyFolderPath+"/"+f)
        
            otherEnvironmentSpace = otherOptions.environmentSpace
            otherParameterSpace   = otherOptions.parameterSpace
        
            environmentSpaceIntersection = set(environmentSpace.keys()).intersection(otherEnvironmentSpace.keys())
            parameterSpaceIntersection   = set(parameterSpace.keys()).intersection(otherParameterSpace.keys())
            if len(set(environmentSpace.keys()))!=len(environmentSpaceIntersection):
                print("ERROR: subprogram failed as environment variables differ from previous experiments found in the output folder.",file=sys.stderr)
                print("environment variables found for CURRENT experiment: " + ", ".join(sorted(environmentSpace.keys())))
                print("environment variables used in PREVIOUS experiment:  " + ", ".join(sorted(otherEnvironmentSpace.keys())))
                sys.exit("Environment variables differ")
            if len(set(parameterSpace.keys()))!=len(parameterSpaceIntersection):
                print("ERROR: subprogram failed as parameters differ from previous experiments found in the output folder.",file=sys.stderr)
                print("parameters found for CURRENT experiment: "+ ", ".join(sorted(parameterSpace.keys())))
                print("parameters used in PREVIOUS experiment:  "+ ", ".join(sorted(otherParameterSpace.keys())))
                sys.exit()

def getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers):
   return scriptsFolderPath + "/" + projectName + "-" + \
          parameterDictHash + "-t"+ranksPerNode+"-c"+cores+"-b"+consumers+".exahype"

def getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,configId,ranks,nodes,run):
   JobInfo = collections.namedtuple("ranksNodesCoreCounts", \
                                    "jobName jobScriptFilePath jobOutputFilePath jobErrorFilePath")
   jobName = projectName + "-" + environmentDictHash + "-" + ungroupedParameterDictHash + \
             "-C"+str(configId)+ "-n" + ranks + "-N" + nodes + "-r"+run
   jobScriptFilePath = scriptsFolderPath + "/" + jobName + ".job"
   jobOutputFilePath = resultsFolderPath + "/" + jobName + ".job_out"
   jobErrorFilePath  = resultsFolderPath + "/" + jobName + ".job_err"
   return JobInfo(jobName,jobScriptFilePath,jobOutputFilePath,jobErrorFilePath) 

def generateScripts():
    """
    Generate spec files and job scripts.
    """
    cpus       = jobs["num_cpus"]
        
    if not os.path.exists(scriptsFolderPath):
        print("create directory "+scriptsFolderPath)
        os.makedirs(scriptsFolderPath)
    
    # spec files
    specFiles=0
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for coreCounts in config.coreCounts:
                cores     = coreCounts.cores
                consumers = coreCounts.consumers
                
                specFilePath = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
                
                specFileBody = renderSpecFile(specFileTemplate,parameterDict,ranksPerNode,cores,consumers)
                with open(specFilePath, "w") as specFile:
                    specFile.write(specFileBody)
                specFiles+=1
    
    print("generated specification files: "+str(specFiles))
    
    # check if required executables exist
    verifyEnvironmentIsCorrect(justWarn=True)
    verifyLogFilterExists(justWarn=True)        
    verifyAllExecutablesExist(justWarn=True)
    
    # generate job scrips
    jobScripts = 0
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)

                    # aggregate the job script body
                    jobScriptBody = ""
                    for coreCount in config.coreCounts:
                        cores     = coreCount.cores
                        consumers = coreCount.consumers
                        for runGrouped in runNumbersGrouped:
                            myRun = run
                            if run=="+":
                                myRun = runGrouped
                            for groupedParameterDict in dictProduct(groupedParameterSpace):
                                parameterDict     = {}
                                parameterDict.update(ungroupedParameterDict)
                                parameterDict.update(groupedParameterDict)
                                parameterDict.pop(None) # ensure we do not hash a dummy None key
                                parameterDictHash = hashDictionary(parameterDict)
                                
                                architecture = parameterDict["architecture"]
                                dimension    = parameterDict["dimension"]
                                ## TODO Continue working here

                                executable   = getExecutableName(environmentDictHash,architecture,dimension,parameterDict)
                                
                                specFilePath   = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
                                                 
                                outputFileName = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                                 "-n" + ranks + "-N" + nodes + "-t"+ranksPerNode+"-c"+cores+"-b"+consumers+"-r"+myRun+".out"
                                outputFilePath = resultsFolderPath + "/" + outputFileName 

                                if "preamble" in jobs:
                                    renderedPreamble = jobs["preamble"].strip("\"")
                                    context = dict(parameterDict)
                                    context["ranks"]=ranks
                                    context["nodes"]=nodes
                                    context["ranksPerNode"]=ranksPerNode
                                    context["coresPerRank"]=cores
                                    context["backgroundTasks"]=consumers
                                    context["run"]=myRun
                                    for key,value in context.items():
                                        renderedPreamble = renderedPreamble.replace("{{"+key+"}}", value)
                                    jobScriptBody += "\n\n" + renderedPreamble + "\n\n"
                                
                                # pipe some information into output file
                                jobScriptBody += "echo \"Timestamp (YYYY/MM/dd:hh:mm:ss): `date +%Y/%m/%d:%H:%M:%S`\" > "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "module list >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "printenv >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "echo \""+jobInfo.jobScriptFilePath+":\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "cat \""+jobInfo.jobScriptFilePath+"\" >> "+outputFilePath+"\n"  
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "echo \""+specFilePath+":\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "cat \""+specFilePath+"\" >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                # pipe environment and parameter dicts into output file
                                jobScriptBody += "echo \"sweep/environment="+json.dumps(environmentDict).replace("\"","\\\"")+"\" >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"sweep/parameters="+json.dumps(parameterDict).replace("\"","\\\"")   +"\" >> "+outputFilePath+"\n"
                                # pipe the commands into the output file
                                runCommand = general["run_command"].replace("\"","")
                                runCommand = runCommand.replace("{{ranks}}",ranks);
                                runCommand = runCommand.replace("{{nodes}}",nodes);
                                runCommand = runCommand.replace("{{ranksPerNode}}",ranksPerNode);
                                runCommand = runCommand.replace("{{cores}}",cores);
                                if "./"==runCommand.strip():
                                    runCommand = runCommand.strip()
                                else:
                                    runCommand += " "
                                jobScriptBody += runCommand+executable+" "+specFilePath+" >> "+outputFilePath+"\n" # no whitespace after runCommand
                                
                                if "likwid" in general:
                                    groups = sweep_options.parseList(general["likwid"])
                                    likwid_pin_arg = int(cores) - 1
                                    if likwid_pin_arg == 0:
                                        likwid_pin_arg = str(likwid_pin_arg)
                                    else:
                                        likwid_pin_arg = "0-"+str(likwid_pin_arg) # e.g. 0-7 pins likwid to cores 0-7 
                                    for group in groups:
                                        likwidCommand = "{runCommand} likwid-perfctr -f -C {likwid_pin_arg} -g {group} {executable} {specFilePath} >> {outputFilePath}.likwid\n"
                                        jobScriptBody += likwidCommand.format(
                                            runCommand=runCommand,
                                            likwid_pin_arg=likwid_pin_arg,
                                            group=group,
                                            executable=executable,
                                            specFilePath=specFilePath,
                                            outputFilePath=outputFilePath
                                        )
                                jobScriptBody += "\n" 
                    
                    # write job file
                    renderedJobScript = renderJobScript(\
                                            jobScriptTemplate,jobScriptBody,jobs,
                                            jobInfo.jobName,jobInfo.jobScriptFilePath,jobInfo.jobOutputFilePath,jobInfo.jobErrorFilePath,
                                            ranks,nodes,ranksPerNode,cores)
                    with open(jobInfo.jobScriptFilePath, "w") as jobScriptFile:
                        jobScriptFile.write(renderedJobScript)
                    
                    jobScripts+=1

    print("generated job scripts: "+str(jobScripts))
                             
def verifyAllJobScriptsExist():
    """
    Verify that all job scripts exist.
    """
    cpus       = jobs["num_cpus"]
    
    if not os.path.exists(scriptsFolderPath):
        print("ERROR: job script folder '"+scriptsFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
    
    allJobScriptsExist = True
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)

                    if not os.path.exists(jobInfo.jobScriptFilePath):
                        allJobScriptsExist = False
                        print("ERROR: job script for " + \
                              "environment="+str(environmentDict)+ \
                              ", (ungrouped)parameters="+str(ungroupedParameterDict) + \
                              ", nodes="+nodes + \
                              ", ranksPerNode="+ranksPerNode + \
                              ", cores="+cores + \
                              ", run="+run + \
                              " does not exist! ('"+jobInfo.jobScriptFilePath+"')",file=sys.stderr)
    if not allJobScriptsExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def verifyAllSpecFilesExist():
    """
    Verify that all ExaHyPE specification files exist.
    """
    cpus       = jobs["num_cpus"]
    
    if not os.path.exists(scriptsFolderPath):
        print("ERROR: job script folder '"+scriptsFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
 
    allSpecFilesExist = True
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for coreCounts in config.coreCounts:
                cores     = coreCounts.cores
                consumers = coreCounts.consumers
                
                specFilePath = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
              
                if not os.path.exists(specFilePath):
                     allSpecFilesExist = False
                     print("ERROR: specification file for \n" + \
                           "parameters="+str(parameterDict) + \
                           ", ranksPerNode="+ranksPerNode + \
                           ", cores="+cores + \
                           " does not exist! ('"+specFilePath+"')",file=sys.stderr)
    
    if not allSpecFilesExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def hashSweep():
    chain = ""
    for config in ranksNodesCoreCounts:
        chain += config.ranks+";"+config.nodes+";"
        for coreCount in config.coreCounts:
            chain += coreCount.cores+";"+coreCount.consumers
    for value in runNumbers:
        chain += str(value)+";"
    for value in runNumbersGrouped:
        chain += str(value)+";"  
        
    for environmentDict in dictProduct(environmentSpace):
        chain += hashDictionary(environmentDict)
    for parameterDict in dictProduct(parameterSpace):
        chain += hashDictionary(parameterDict)
        
    return hashlib.md5(chain.encode()).hexdigest()[0:8]

def extractJobId(processOutput):
    jobId = "unknown"
    lines = processOutput.split("\n")
    for line in lines:
        # SLURM
        # hamilton: 'Submitted batch job 67586'
        # coolmuc:  'Submitted batch job 67586 on cluster mpp3'
        # LOAD-LEVELER:
        # supermuc: 'llsubmit: The job "srv24ib.840220" has been submitted.'
        if "Submitted batch job " in line:
            jobId = line.strip().split(" ")[3]
        if "llsubmit: The job " in line:
            jobId = line.strip().split("\"")[1]
    return jobId

def submitJobs(submitAllJobs=False):
    """
    Submit all jobs spanned by the options.
    """
    jobSubmissionTool    = general["job_submission"].replace("\"","")

    jobIdsPath       = historyFolderPath + "/" + hashSweep() + ".submitted" # TODO rename file
    acceptedJobsPath = historyFolderPath + "/" + hashSweep() + ".accepted-jobs"
    
    cpus = jobs["num_cpus"]
    
    # verify everything is fine
    verifyEnvironmentIsCorrect()
    verifyLogFilterExists(justWarn=True)        
    verifyAllExecutablesExist()
    verifyAllJobScriptsExist()
    verifyAllSpecFilesExist()
    
    if not os.path.exists(resultsFolderPath):
        print("create directory "+resultsFolderPath)
        os.makedirs(resultsFolderPath)

    jobIds       = []
    if os.path.exists(jobIdsPath):
        with open(jobIdsPath, "r") as file:
            jobIds.extend(json.loads(file.read()))

    acceptedJobs = []
    if submitAllJobs and os.path.exists(acceptedJobsPath):
        with open(acceptedJobsPath, "r") as file:
            acceptedJobs.extend(json.loads(file.read()))

    # loop over job scrips
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)
                            
                    command=jobSubmissionTool + " " + jobInfo.jobScriptFilePath
                    if command not in acceptedJobs:
                        print(command)
                        process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, err) = process.communicate()
                        process.wait()
                        if err:
                            print(err)
                        else:
                            acceptedJobs.append(command)
                            jobIds.append(extractJobId(output.decode("UTF_8")))
    
    if not os.path.exists(historyFolderPath):
        print("create directory "+historyFolderPath)
        os.makedirs(historyFolderPath)
    
    with open(jobIdsPath, "w") as file:
        file.write(json.dumps(jobIds))
    with open(acceptedJobsPath, "w") as file:
        file.write(json.dumps(acceptedJobs))
    
    print("submitted "+str(len(jobIds))+" jobs")
    print("job ids are memorised in: "+jobIdsPath)
    print("accepted job scripts are stored in: "+acceptedJobsPath)
    command="cp "+optionsFile+" "+jobIdsPath.replace(".submitted",".ini")
    print(command)
    process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    (output, err) = process.communicate()
    process.wait()

def cancelJobs():
    """
    Cancel submitted jobs.
    """
    jobCancellationTool = general["job_cancellation"].replace("\"","")
    
    jobIdsPath       = historyFolderPath + "/" + hashSweep() + ".submitted" # TODO rename file
    acceptedJobsPath = historyFolderPath + "/" + hashSweep() + ".accepted-jobs"
    
    jobIds = None
    try:
        with open(jobIdsPath, "r") as jobIdsFile:
            jobIds = json.loads(jobIdsFile.read())
    except IOError:
        print("ERROR: couldn't find any submitted jobs for current sweep ('"+jobIdsPath+"').")
        sys.exit()

    for jobId in jobIds:
        command = jobCancellationTool + " " + jobId
        print(command)
        subprocess.call(command,shell=True)
    print("cancelled "+str(len(jobIds))+" jobs")    

    os.remove(jobIdsPath)
    os.remove(acceptedJobsPath)

if __name__ == "__main__":
    import sys,os
    import subprocess
    import itertools
    import hashlib
    import json
    import re
    import math
    import collections   
 
    import sweep_analysis
    import sweep_options
    
    subprograms = [\
"build","buildMissing","buildLocally","link","scripts","submit","submitAll","cancel","parseAdapters",\
"parseTotalTimes","parseTimeStepTimes","parseMetrics","parseJobStatistics",\
"cleanBuild", "cleanScripts","cleanResults","cleanHistory","cleanAll",\
"iniTemplate","jobTemplate"]
    
    if haveToPrintHelpMessage(sys.argv):
        info = \
"""sweep.py:

1) run:

./sweep.py myoptions.ini <subprogram>

available subprograms:

* build                 - build all executables
* buildMissing          - build only missing executables
* buildLocally          - rebuild only the local application folder (no make clean)
* link                  - link runtime dependencies into the build folder
* scripts               - generate specification files and job scripts
* submit                - submit only the generated job scripts, which have not been submitted yet or have been rejected by the scheduler
* submitAll             - submit all jobs
* cancel                - cancel the submitted jobs
* parseAdapters         - read the job output and parse adapter times
* parseTotalTimes       - read the adapter times table and:
                          Per configuration, calculate the total simulation time 
                          and minimise over all runs.
* parseTimeStepTimes    - read the adapter times table and:
                          Per configuration, calculate the time spent per time step
                          and minimise over all runs.
* parseMetrics          - read the job output and parse likwid metrics
* parseJobStatistics    - read background job processing stastistics. 
                          Requires that executables are built with "-DTBB_USE_THREADING_TOOLS=1".
* cleanAll              - remove the whole sweep benchmark suite
* cleanBuild            - remove the build subfolder
* cleanScripts          - remove the scripts subfolder
* cleanResults          - remove the results subfolder
* cleanHistory          - clean the submission history
* iniTemplate          - print out template for ini file
* jobTemplate <machine> - print out job template file for supercomputer. Supported: "supermuc", "supermuc-skx", "hamilton"

2) typical workflow:

./sweep.py iniTemplate > myTemplateFile.ini
./sweep.py jobTemplate supermuc > myJobTemplate.ini

# edit both files and create specification file template

./sweep.py myoptions.ini build
./sweep.py myoptions.ini scripts
./sweep.py myoptions.ini submit

(after jobs have finished)

./sweep.py myoptions.ini parseAdapters [--compress]
./sweep.py myoptions.ini parseTotalTimes
./sweep.py myoptions.ini parseTimeStepTimes

note: with --compress you remove table columns containing the
same value in every row.

(if likwid measurements have been performed)

./sweep.py myoptions.ini parseMetrics [--compress]


3) general options file (*.ini) structure

the options file must contain the following sections:

* [general]     - general settings such as the ExaHyPE-Engine path, the number of make threads, how to submit and cancel jobs, and so on
* [jobs]        - settings for the jobs to be run, e.g. the number of nodes, ranksPerNode, cores, and so on
* [environment] - environments for the executables to be built, e.g. the MODE, SHAREDMEM, DISTRIBUTEDMEM environment variables can be specified here

It must further contain at least one of the following sections:

* [parameters]         - The specification file parameters. Most of them are optional but some are required such as 'dimension', 'order', and so on. 
* [parameters_grouped] - More specificiation file parameters. But the corresponding runs are placed all
                         on the same job script.

"""
        print(info) # correctly indented
        sys.exit()
    
    optionsFile = parseArgument(sys.argv,1)
    
    # print out templates
    if optionsFile == "iniTemplate":
        sweep_options.printOptionsFileTemplate()
        sys.exit()
    elif optionsFile == "jobTemplate":
        machine = parseArgument(sys.argv,2)
        sweep_options.printJobTemplate(machine)
        sys.exit()

    # run subprograms
    subprogram  = parseArgument(sys.argv,2)
    compressTable = parseArgument(sys.argv,3)=="--compress"

    options = sweep_options.parseOptionsFile(optionsFile)
    
    general                   = options.general
    jobs                      = options.jobs
    environmentSpace          = options.environmentSpace
    parameterSpace            = options.parameterSpace
    ungroupedParameterSpace   = options.ungroupedParameterSpace
    groupedParameterSpace     = options.groupedParameterSpace
    compileTimeParameterSpace = options.compileTimeParameterSpace
 
    exahypeRoot      = options.exahypeRoot
    outputPath       = options.outputPath
    projectPath      = options.projectPath
    projectName      = options.projectName
   
    buildFolder       = options.buildFolder
    buildFolderPath   = options.buildFolderPath
    scriptsFolderPath = options.scriptsFolderPath
    resultsFolderPath = options.resultsFolderPath
    historyFolderPath = options.historyFolderPath
    
    jobClass             = options.jobClass
    islands              = options.islands
    ranksNodesCoreCounts = options.ranksNodesCoreCounts
    runNumbers           = options.runNumbers
    runNumbersGrouped    = options.runNumbersGrouped   

    createdFirstSpecFile  = False
    createdFirstJobScript = False
 
    verifySweepAgreesWithHistoricalExperiments()
    
    specFileTemplatePath = exahypeRoot+"/"+general["spec_template"]
    specFileTemplate     = None
    try:
        with open(specFileTemplatePath, "r") as specFileTemplateFile:
            specFileTemplate=specFileTemplateFile.read()
    except IOError:
        print("ERROR: couldn\'t open specification file template file: "+specFileTemplatePath,file=sys.stderr)
        sys.exit("Specification file missing")
        
    jobScriptTemplatePath = exahypeRoot+"/"+general["job_template"]    
    jobScriptTemplate = None
    try:
        with open(jobScriptTemplatePath, "r") as jobScriptTemplateFile:
            jobScriptTemplate=jobScriptTemplateFile.read()
    except IOError:
        print("ERROR: couldn\'t open job script template file: "+jobScriptTemplatePath,file=sys.stderr)
        sys.exit("Job Template missing")
    
    # TODO move into options?
    verifyAllRequiredParametersAreGiven(specFileTemplate)
 
    # select subprogram
    if subprogram == "cleanAll":
        clean()
    elif subprogram == "cleanBuild":
        clean("build")
    elif subprogram == "cleanScripts":
        clean("scripts")
    elif subprogram == "cleanResults":
        clean("results")
    elif subprogram == "cleanHistory":
        clean("history")
    elif subprogram == "build":
        build(buildOnlyMissing=False, skipMakeClean=False)
    elif subprogram == "buildMissing":
        build(buildOnlyMissing=True, skipMakeClean=False)
    elif subprogram == "buildLocally":
        build(buildOnlyMissing=False, skipMakeClean=True)
    elif subprogram == "link":
        unlink()
        link()
    elif subprogram == "scripts":
        generateScripts()
    elif subprogram == "submit":
        submitJobs(submitAllJobs=False)
    elif subprogram == "submitJob":
        submitJobs(submitAllJobs=True)
    elif subprogram == "cancel":
        cancelJobs()
    elif subprogram == "parseAdapters":
        sweep_analysis.parseAdapterTimes(resultsFolderPath,projectName,compressTable)
    elif subprogram == "parseTotalTimes":
        sweep_analysis.parseSummedTimes(resultsFolderPath,projectName)
    elif subprogram == "parseTimeStepTimes":
        sweep_analysis.parseSummedTimes(resultsFolderPath,projectName,timePerTimeStep=True)
    elif subprogram == "parseMetrics":
        sweep_analysis.parseMetrics(resultsFolderPath,projectName,compressTable)
    elif subprogram == "parseJobStatistics":
        sweep_analysis.parseJobStatistics(resultsFolderPath,projectName,compressTable)
