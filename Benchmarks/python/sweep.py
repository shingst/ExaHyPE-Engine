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
    parameterSpace = {}
    if "parameters" in config and len(config["parameters"].keys()):
        for key, value in config["parameters"].items():
            parameterSpace[key] = parseList(value)
            
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
    else:
        print("ERROR: Section 'parameters' is missing or empty! (Must contain at least 'dimension' and 'order'.)",file=sys.stderr)
        sys.exit()
    
    return parameterSpace

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
    
    return hashlib.md5(chain.encode()).hexdigest()

def clean(subFolder=""):
    """
    Clean the complete output folder or just a subfolder
    if specified.
    """
    folder = exahypeRoot+"/"+outputPath+"/"+subFolder
    print("rm -r "+folder)
    subprocess.call("rm -r "+folder, shell=True)

def renderSpecFile(templateBody,parameterDict,tasks,cores):
    renderedFile = templateBody
    
    context = dict(parameterDict)
    
    consistent = True
    # verify mandatory options file parameters can be found in template file
    keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",templateBody)]
    for key in parameterDict:
        if key not in keysInTemplate:
            consistent = False
            print("ERROR: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr) 
    # optional parameters
    context["tasks"] = tasks
    context["cores"] = cores
    # verify template parameters are defined in options file
    for key in keysInTemplate:
        if key not in context:
            consistent = False
            print("ERROR: specification file template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
    if not consistent:
        print("ERROR: subprogram aborted as specification file template and sweep options file are inconsistent.",file=sys.stderr)
        sys.exit()
    
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def build(buildOnlyMissing=False):
    """
    Build the executables.
    
    Args:
    buildOnlyMissing(bool):
       Build only missing executables.
    """
    print("currently loaded modules:")
    subprocess.call("module list",shell=True)
    print("")
    print("ExaHyPE build environment:")
    exahypeEnv = ["COMPILER", "MODE", "SHAREDMEM", "DISTRIBUTEDMEM", "EXAHYPE_CC", "PROJECT_C_FLAGS", "PROJECT_L_FLAGS", "COMPILER_CFLAGS", "COMPILER_LFLAGS", "FCOMPILER_CFLAGS", "FCOMPILER_LFLAGS"]
    for variable in exahypeEnv:
        if variable in os.environ:
            print(variable+"="+os.environ[variable])
    print("")
    
    templateFileName = general["spec_template"]
    
    templateBody = None
    try:
        with open(exahypeRoot+"/"+templateFileName, "r") as templateFile:
            templateBody=templateFile.read()
    except IOError:
        print("ERROR: couldn\'t open template file: "+templateFileName,file=sys.stderr)
        sys.exit()
    
    if not os.path.exists(buildFolderPath):
        print("create directory "+buildFolderPath)
        os.makedirs(buildFolderPath)
        
    architectures = parameterSpace["architecture"]
    optimisations = parameterSpace["optimisation"]
    dimensions    = parameterSpace["dimension"]
    orders        = parameterSpace["order"]
    buildParameterDict = list(dictProduct(parameterSpace))[0]
        
    firstIteration = True
    executables=0
    for environmentDict in dictProduct(environmentSpace):
        for key,value in environmentDict.items():
            os.environ[key]=value
        environmentDictHash = hashDictionary(environmentDict)   
        
        for architecture in architectures:
            for optimisation in optimisations:
                for dimension in dimensions:
                    if not firstIteration:
                        command = "make clean"
                        print(command)
                        process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        (output, err) = process.communicate()
                        process.wait()
                    for order in orders:
                        oldExecutable = exahypeRoot + "/" + projectPath+"/ExaHyPE-"+projectName
                        executable = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+\
                                     architecture+"-"+optimisation+"-d" + dimension + "-p" + order
                            
                        if not os.path.exists(executable) or not buildOnlyMissing:
                            buildParameterDict["optimisation"]=optimisation
                            buildParameterDict["architecture"]=architecture
                            buildParameterDict["dimension"]   =dimension
                            buildParameterDict["order"]       =order
                                
                            buildSpecFileBody = renderSpecFile(templateBody,buildParameterDict,"1","1")
                                
                            buildSpecFilePath = outputPath+"/"+buildFolder+"/"+projectName+"-"+architecture+"-"+optimisation+"-d"+dimension+"-p"+order+".exahype"
                                
                            with open(exahypeRoot + "/" + buildSpecFilePath, "w") as buildSpecificationFile:
                                buildSpecificationFile.write(buildSpecFileBody)
                                
                            print("building executable for " + \
                                  "environment="+str(environmentDict) + \
                                  ", architecture="+architecture + \
                                  ", optimisation="+optimisation + \
                                  ", dimension="+dimension + \
                                  ", order="+order,\
                                  file=sys.stderr)
                            # run toolkit
                            toolkitCommand = "(cd "+exahypeRoot+" && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive "+buildSpecFilePath+")"
                            print(toolkitCommand,end="",flush=True)
                            process = subprocess.Popen([toolkitCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            (output, toolkitErr) = process.communicate()
                            process.wait()
                            if "setup build environment ... ok" in str(output):
                                print(" [OK]")
                            else:
                                print(" [FAILED]")
                                print("toolkit output=\n"+output.decode('UTF-8'),file=sys.stderr)
                                print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                                sys.exit()

                            # clean application folder only
                            command = "rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels"
                            print(command)
                            process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            process.communicate()
                            process.wait()
                            
                            if firstIteration:
                                command = "make clean"
                                print(command)
                                process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                (output, err) = process.communicate()
                                process.wait()
                                firstIteration = False

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
                                print(" [FAILED]")
                                print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                                sys.exit()

                            moveCommand   = "mv "+oldExecutable+" "+executable
                            print(moveCommand)
                            subprocess.call(moveCommand,shell=True)
                            print("--------------------------------------------------------------------------------")
                            print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                            print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                            print("--------------------------------------------------------------------------------")
                            executables+=1
                        else:
                            print("skipped building of '"+executable+"' as it already exists.")

    print("built executables: "+str(executables))


def renderJobScript(templateBody,environmentDict,parameterDict,jobs,
                    jobName,jobFilePath,outputFileName,errorFileName,appName,specFilePath,
                    nodes,tasks,cores,run):
    """
    Render a job script.
    """
    renderedFile = templateBody
    
    context = {}
    # mandatory
    context["nodes"] = nodes
    context["tasks"] = tasks
    context["cores"] = cores
    context["job_name"]    = jobName
    context["output_file"] = outputFileName
    context["error_file"]  = errorFileName
    
    context["environment"] = json.dumps(environmentDict).replace("\"","\\\"")
    context["parameters"]  = json.dumps(parameterDict).replace("\"","\\\"")
    
    context["job_file"]  = jobFilePath
    context["app"]       = appName
    context["spec_file"] = specFilePath
    
    consistent = True
    # verify all mandatory(!) sweep options are defined in template
    keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",templateBody)]
    for key in context:
        if key not in keysInTemplate:
            consistent = False
            print("ERROR: parameter '{{"+key+"}}' not found in job script template!",file=sys.stderr)
    
    # put optional sweep options in context
    context["mail"]  = jobs["mail"]
    context["time"]  = jobs["time"]
    context["ranks"] = str(int(nodes)*int(tasks))
    
    # now verify template parameters are defined in options file
    for key in keysInTemplate:
        if key not in context:
            consistent = False
            print("ERROR: job script template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
    if not consistent:
        print("ERROR: subprogram aborted as job script template and sweep options file are inconsistent.",file=sys.stderr)
        sys.exit()
    
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def verifyAllExecutablesExist(justWarn=False):
    """
    Verify that all executables exist.
    """
    architectures = parameterSpace["architecture"]
    optimisations = parameterSpace["optimisation"]
    dimensions    = parameterSpace["dimension"]
    orders        = parameterSpace["order"]
    
    messageType = "ERROR"
    if justWarn:
      messageType = "WARNING"
    
    if not justWarn and not os.path.exists(buildFolderPath):
        print("ERROR: build folder '"+buildFolderPath+"' doesn't exist! Please run subprogram 'build' beforehand.",file=sys.stderr)
        sys.exit()
    
    allExecutablesExist = True
    for environmentDict in dictProduct(environmentSpace):
        environmentDictHash = hashDictionary(environmentDict)
        for architecture in architectures:
            for optimisation in optimisations:
                for dimension in dimensions:
                    for order in orders:
                        executable = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+\
                                     architecture+"-"+optimisation+"-d" + dimension + "-p" + order
                
                        if not os.path.exists(executable):
                            allExecutablesExist = False
                            print(messageType+ ": application for " + \
                                  "environment="+str(environmentDict) + \
                                  ", dimension="+dimension + \
                                  ", order="+order + \
                                  " does not exist! ('"+executable+"')",file=sys.stderr)
            
    if not justWarn and not allExecutablesExist:
        print("ERROR: subprogram failed as not all executables exist. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the 'build' subprogram.",file=sys.stderr)
        sys.exit()

def verifySweepAgreesWithHistoricExperiments():
    """
    If there are any previous experiments ensure that the sweep 
    parameter spaces contain the same axes.
    """
    if os.path.exists(historyFolderPath):
        previousSweeps = [f for f in os.listdir(historyFolderPath) if f.endswith(".ini")]
        print(previousSweeps)
    
        for f in previousSweeps:
            otherConfig = configparser.ConfigParser()
            otherConfig.optionxform=str
            otherConfig.read(historyFolderPath + "/" + f)
        
            otherEnvironmentSpace = parseEnvironment(otherConfig)
            otherParameterSpace   = parseParameters(otherConfig)
        
            environmentSpaceIntersection = set(environmentSpace.keys()).intersection(otherEnvironmentSpace.keys())
            parameterSpaceIntersection   = set(parameterSpace.keys()).intersection(otherParameterSpace.keys())
            if len(set(environmentSpace.keys()))!=len(environmentSpaceIntersection):
                print("ERROR: subprogram failed as environment variables differ from previous experiments found in the output folder.",file=sys.stderr)
                print("environment variables found for CURRENT experiment: " + ", ".join(sorted(environmentSpace.keys())))
                print("environment variables used in PREVIOUS experiment:  " + ", ".join(sorted(otherEnvironmentSpace.keys())))
                sys.exit()
            if len(set(parameterSpace.keys()))!=len(parameterSpaceIntersection):
                print("ERROR: subprogram failed as parameters differ from previous experiments found in the output folder.",file=sys.stderr)
                print("parameters found for CURRENT experiment: "+ ", ".join(sorted(parameterSpace.keys())))
                print("parameters used in PREVIOUS experiment:  "+ ", ".join(sorted(otherParameterSpace.keys())))
                sys.exit()

def generateScripts():
    """
    Generate spec files and job scripts.
    """
    cpus       = jobs["num_cpus"]
    
    specFileTemplatePath  = exahypeRoot+"/"+general["spec_template"]
    jobScriptTemplatePath = exahypeRoot+"/"+general["job_template"]
    
    specFileTemplate  = None
    try:
        with open(specFileTemplatePath, "r") as templateFile:
            specFileTemplate=templateFile.read()
    except IOError: 
        print("ERROR: couldn\'t open template file: "+specFileTemplatePath,file=sys.stderr)
        sys.exit()
        
    jobScriptTemplate = None
    try:
        with open(jobScriptTemplatePath, "r") as templateFile:
            jobScriptTemplate=templateFile.read()
    except IOError:
        print("ERROR: couldn\'t open template file: "+jobScriptTemplatePath,file=sys.stderr)
        sys.exit()
        
    if not os.path.exists(scriptsFolderPath):
        print("create directory "+scriptsFolderPath)
        os.makedirs(scriptsFolderPath)
    
    # spec files
    specFiles=0
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for tasks in taskCounts:
            for parsedCores in coreCounts:
              cores = parsedCores
              if parsedCores=="auto":
                 cores=str(int(int(cpus) / int(tasks)))
              specFileBody = renderSpecFile(specFileTemplate,parameterDict,tasks,cores)
              
              specFilePath = scriptsFolderPath + "/" + projectName + "-" + parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
              
              with open(specFilePath, "w") as specFile:
                  specFile.write(specFileBody)
              specFiles+=1
    
    print("generated specification files: "+str(specFiles))
    
    # check if required executables exist
    verifyAllExecutablesExist(True)
    
    # generate job scrips
    jobScripts = 0
    for run in runNumbers:
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for parsedCores in coreCounts:
                    cores = parsedCores
                    if parsedCores=="auto":
                        cores=str(int(int(cpus) / int(tasks)))
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            architecture = parameterDict["architecture"]
                            optimisation = parameterDict["optimisation"]
                            dimension    = parameterDict["dimension"]
                            order        = parameterDict["order"]
                            
                            executable   = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+\
                                           architecture +"-"+ optimisation +"-d" + dimension + "-p" + order
                            specFilePath = scriptsFolderPath + "/" + projectName + "-" + \
                                           parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+run
                            jobFilePrefix  = scriptsFolderPath + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            outputFileName = resultsFolderPath + "/" + jobName + ".out"
                            errorFileName  = resultsFolderPath + "/" + jobName + ".err"
                            
                            jobScriptBody = renderJobScript(jobScriptTemplate,environmentDict,parameterDict,jobs,
                                                            jobName,jobFilePath,outputFileName,errorFileName,executable,specFilePath,
                                                            nodes,tasks,cores,run)
                            with open(jobFilePath, "w") as jobFile:
                                jobFile.write(jobScriptBody)
                            
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
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for parsedCores in coreCounts:
                    cores = parsedCores
                    if parsedCores=="auto":
                        cores=str(int(int(cpus) / int(tasks)))
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+run
                            jobFilePrefix  = scriptsFolderPath + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            
                            if not os.path.exists(jobFilePath):
                                allJobScriptsExist = False
                                print("ERROR: job script for " + \
                                      "environment="+str(environmentDict)+ \
                                      ", parameters="+str(parameterDict) + \
                                      ", nodes="+nodes + \
                                      ", tasks="+tasks + \
                                      ", cores="+cores + \
                                      ", run="+run + \
                                      " does not exist! ('"+jobFilePath+"')",file=sys.stderr)
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
        
        for tasks in taskCounts:
            for parsedCores in coreCounts:
                cores = parsedCores
                if parsedCores=="auto":
                    cores=str(int(int(cpus) / int(tasks)))
                
                specFilePath = scriptsFolderPath + "/" + projectName + "-" + parameterDictHash + "-t"+tasks+"-c"+cores+".exahype"
              
                if not os.path.exists(specFilePath):
                     allSpecFilesExist = False
                     print("ERROR: specification file for \n" + \
                           "parameters="+str(parameterDict) + \
                           ", tasks="+tasks + \
                           ", cores="+cores + \
                           " does not exist! ('"+specFilePath+"')",file=sys.stderr)
    
    if not allSpecFilesExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def hashSweep():
    chain = ""
    for value in nodeCounts:
        chain += value+";"
    for value in taskCounts:
        chain += value+";"
    for value in coreCounts:
        chain += value+";"
    for value in runNumbers:
        chain += value+";"
    
    for environmentDict in dictProduct(environmentSpace):
        chain += hashDictionary(environmentDict)
    for parameterDict in dictProduct(parameterSpace):
        chain += hashDictionary(parameterDict)
        
    return hashlib.md5(chain.encode()).hexdigest()

def extractJobId(processOutput):
    jobId = "unknown"
    lines = processOutput.split("\n")
    for line in lines:
        # SLURM
        # hamilton: "Submitted batch job 67586"
        # coolmuc:  "Submitted batch job 67586 on cluster mpp3"
        if "Submitted batch job " in line:
            jobId = line.strip().split(" ")[3]
    return jobId

def submitJobs():
    """
    Submit all jobs spanned by the options.
    """
    jobSubmissionTool    = general["job_submission"]
    
    cpus = jobs["num_cpus"]
    
    # verify everything is fine
    verifyAllExecutablesExist()
    verifyAllJobScriptsExist()
    verifyAllSpecFilesExist()
    
    if not os.path.exists(resultsFolderPath):
        print("create directory "+resultsFolderPath)
        os.makedirs(resultsFolderPath)
    
    # loop over job scrips
    jobIds = []
    for run in runNumbers:
        for nodes in nodeCounts:
            for tasks in taskCounts:
                for parsedCores in coreCounts:
                    cores = parsedCores
                    if parsedCores=="auto":
                        cores=str(int(int(cpus) / int(tasks)))
                    for environmentDict in dictProduct(environmentSpace):
                        environmentDictHash = hashDictionary(environmentDict)
                        
                        for parameterDict in dictProduct(parameterSpace):
                            parameterDictHash = hashDictionary(parameterDict)
                            
                            dimension = parameterDict["dimension"]
                            order     = parameterDict["order"]
                            
                            jobName        = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                             "-n" + nodes + "-t"+tasks+"-c"+cores+"-r"+run
                            jobFilePrefix  = scriptsFolderPath + "/" + jobName
                            jobFilePath    = jobFilePrefix + ".job"
                            
                            command=jobSubmissionTool + " " + jobFilePath
                            print(command)
                            process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                            (output, err) = process.communicate()
                            process.wait()
                            jobIds.append(extractJobId(output.decode("UTF_8")))
    
    if not os.path.exists(historyFolderPath):
        print("create directory "+historyFolderPath)
        os.makedirs(historyFolderPath)
    
    submittedJobsPath = historyFolderPath + "/" + hashSweep() + ".submitted"
    
    with open(submittedJobsPath, "w") as submittedJobsFile:
        submittedJobsFile.write(json.dumps(jobIds))
    
    print("submitted "+str(len(jobIds))+" jobs")
    print("job ids are memorised in: "+submittedJobsPath)
    command="cp "+configFile+" "+submittedJobsPath.replace(".submitted",".ini")
    print(command)
    process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    (output, err) = process.communicate()
    process.wait()

def cancelJobs():
    """
    Cancel submitted jobs.
    """
    jobCancellationTool    = general["job_cancellation"]
    
    submittedJobsPath = historyFolderPath + "/" + hashSweep() + ".submitted"
    
    jobIds = None
    try:
        with open(submittedJobsPath, "r") as submittedJobsFile:
            jobIds = json.loads(submittedJobsFile.read())
    except IOError:
        print("ERROR: couldn't find any submitted jobs for current sweep ('"+submittedJobsPath+"').")
        sys.exit()

    for jobId in jobIds:
        command = jobCancellationTool + " " + jobId
        print(command)
        subprocess.call(command,shell=True)
    print("cancelled "+str(len(jobIds))+" jobs")    

    command = "rm "+submittedJobsPath
    command = "rm "+submittedJobsPath.replace(".submitted",".ini")
    print(command)
    subprocess.call(command,shell=True)

def parseResultFile(filePath):
    '''
    Reads a single sweep job output file and parses the user time spent within each adapter.
    
    Args:
       filePath (str):
          Path to the sweep job output file.

    Returns:
       A dict holding for each of the found adapters a nested dict that holds the following key-value pairs:
          * 'n'       : (int)    Number of times this adapter was used.
          * 'cputime' : (float) Total CPU time spent within the adapter.
          * 'usertime': (float) Total user time spent within the adapter.
       
       The dict further holds the following dictionaries:
          * 'environment':(dictionary(str,str)) Total user time spent within the adapter.
          * 'parameters' :(dictionary(str,str)) Total user time spent within the adapter.
    '''
    environmentDict = {}
    parameterDict   = {}
    
    adapters = {}
    cputimeIndex  = 3
    usertimeIndex = 5
    
    try:
        fileHandle=codecs.open(filePath,'r','UTF_8')
        for line in fileHandle:
            if line.startswith("sweep/environment"):
                value = line.split('=')[-1]
                environmentDict=json.loads(value)
            if line.startswith("sweep/parameters"):
                value = line.split('=')[-1]
                parameterDict=json.loads(value)
            anchor = '|'
            header = '||'
            if anchor in line and header not in line:
                segments = line.split('|')
                adapter = segments[1].strip();
                adapters[adapter]             = {}
                adapters[adapter]['iterations']     = segments[2].strip()
                adapters[adapter]['total_cputime']  = segments[cputimeIndex ].strip()
                adapters[adapter]['total_usertime'] = segments[usertimeIndex].strip()
    except IOError as err:
        print ("ERROR: could not parse adapter times for file "+filePath+"! Reason: "+str(err))
    return environmentDict,parameterDict,adapters

def getAdapterTimesSortingKey(row):
    keyTuple = ()
    keys = row[:-3]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseAdapterTimes():
    """
    Loop over all ".out" files in the results section and create a table.
    """
    tablePath         = resultsFolderPath+"/"+projectName+'.csv'
    try:
        with open(tablePath, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=';')
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out")]

            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-t1-c1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-t([0-9]+)-c([0-9]+)-r([0-9]+).out$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                nodes               = match.group(4)
                tasks               = match.group(5)
                cores               = match.group(6)
                run                 = match.group(7)
                
                environmentDict,parameterDict,adapters = parseResultFile(resultsFolderPath + "/" + fileName)
                if len(adapters):
                    # write header
                    if firstFile:
                        header = []
                        header += sorted(environmentDict)
                        for parameter in knownParameters:
                            header.append(parameter)
                        for parameter in sorted(parameterDict):
                            if parameter not in knownParameters:
                                header.append(parameter)
                        header.append("nodes")
                        header.append("tasks")
                        header.append("cores")
                        header.append("adapter")
                        header.append("run")
                        header.append("iterations")
                        header.append("total_cputime")
                        header.append("total_usertime")
                        csvwriter.writerow(header)
                        firstFile=False
                    print(resultsFolderPath+"/"+fileName)
                    
                    # write rows
                    for adapter in adapters:
                        row=[]
                        for key in sorted(environmentDict):
                            row.append(environmentDict[key])
                        for key in knownParameters:
                            row.append(parameterDict[key])
                        for key in sorted(parameterDict):
                            if key not in knownParameters:
                                row.append(parameterDict[key])
                        row.append(nodes)
                        row.append(tasks)
                        row.append(cores)
                        row.append(adapter)
                        row.append(run)
                        row.append(adapters[adapter]["iterations"])
                        row.append(adapters[adapter]["total_cputime"])
                        row.append(adapters[adapter]["total_usertime"])
                        csvwriter.writerow(row)
        success = not firstFile
        if success:
          # reopen the file and sort it
          tableFile   = open(tablePath, 'r')
          header      = next(tableFile)
          header      = header.strip()
          reader      = csv.reader(tableFile,delimiter=',')
          
          sortedData = sorted(reader,key=getAdapterTimesSortingKey)
          tableFile.close()
          
          with open(tablePath, 'w') as sortedTableFile:
              writer = csv.writer(sortedTableFile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
              writer.writerow(header.split(','))
              writer.writerows(sortedData)
          print("created table:")
          print(tablePath) 

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: "<<str(err))


def parseLikwidMetrics(filePath,singlecore=False):
    """
    Reads a single Peano output file and parses likwid performance metrics.
    
    Args:
       filePath (str):
          Path to the Peano output file.
       metrics (str[][]):
          A list of metrics the we want to read out.
       counters (str[][]):
          A list of counters the we want to read out.
       singlecore (bool):
          Specifies if the run was a singlecore run.
    
    Returns:
       A dict holding for each of the found metrics a nested dict that holds the following key-value pairs:
          * 'Sum' 
          * 'Avg' 
          * 'Min' 
          * 'Max' 
    """
    columns    = [ "Sum","Min","Max","Avg" ]
    
    environmentDict = {}
    parameterDict   = {}
    
    result  = {}
    for metric in metrics:
        result[metric[0]] =  {}
        result[metric[0]][metric[1]] = -1.0
    for counter in counters:
        result[counter[0]] =  {}
        result[counter[0]][counter[1]] = -1.0

    try:
        fileHandle=open(filePath)
        
        for line in fileHandle:
            if line.startswith("sweep/environment"):
                value = line.split('=')[-1]
                environmentDict=json.loads(value)
            if line.startswith("sweep/parameters"):
                value = line.split('=')[-1]
                parameterDict=json.loads(value)
            
            for metric in metrics: 
                if singlecore:
                    if metric[0] in line:
                        segments = line.split('|')
                        
                        #    |     Runtime (RDTSC) [s]    |    6.5219    |
                        value  = float(segments[2].strip());
                        values = {}                         
                        values["Sum"] = value
                        values["Min"] = value
                        values["Max"] = value
                        values["Avg"] = value
                        result[metric[0]][metric[1]]=values[metric[1]]                        
                else:
                    if metric[0]+" STAT" in line:
                        segments = line.split('|')
                        #   |  Runtime (RDTSC) [s] STAT |   27.4632  |   1.1443  |   1.1443  |   1.1443  |
                        values = {}                                                 
                        values["Sum"] = float(segments[2].strip());
                        values["Min"] = float(segments[3].strip());
                        values["Max"] = float(segments[4].strip());
                        values["Avg"] = float(segments[5].strip());
                        result[metric[0]][metric[1]]=values[metric[1]]
                        
            for counter in counters: 
                if singlecore:
                    if counter[0] in line:
                        segments = line.split('|')
                        #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE   |   PMC1  |  623010105225  | ...
                        value  = float(segments[3].strip());
                        values = {}                         
                        values["Sum"] = value
                        values["Min"] = value
                        values["Max"] = value
                        values["Avg"] = value
                        result[counter[0]][counter[1]]=values[metric[1]]                        
                else:
                    if counter[0]+" STAT" in line:
                        segments = line.split('|')
                        #    |    FP_ARITH_INST_RETIRED_SCALAR_DOUBLE STAT   |   PMC1  |  623010105225  | ...
                        values = {}                                                 
                        values["Sum"] = float(segments[3].strip());
                        values["Min"] = float(segments[4].strip());
                        values["Max"] = float(segments[5].strip());
                        values["Avg"] = float(segments[6].strip());
                        result[counter[0]][counter[1]]=values[counter[1]]
    except IOError as err:
        print ("ERROR: could not parse likwid metrics for file "+filePath+"! Reason: "+str(err))
    return environmentDict,parameterDict,result

def getLikwidMetricsSortingKey(row):
    keyTuple = ()
    keys = row[:-(len(metrics)+len(counters))]
    for key in keys:
      try:
          keyTuple += (float(key),)
      except ValueError:
          keyTuple += (key,)
    return keyTuple

def parseMetrics():
    """
    Loop over all ".out.likwid" files in the results section and create a table.
    """
    tablePath         = resultsFolderPath+"/"+projectName+'-likwid.csv'
    try:
        with open(tablePath, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=';')
            files = [f for f in os.listdir(resultsFolderPath) if f.endswith(".out.likwid")]

            print("processed files:")
            firstFile = True
            for fileName in files:
                # example: Euler-088f94514ee5a8f92076289bf648454e-26b5e7ccb0354b843aad07aa61fd110d-n1-t1-c1-r1.out
                match = re.search('^(.+)-(.+)-(.+)-n([0-9]+)-t([0-9]+)-c([0-9]+)-r([0-9]+).out.likwid$',fileName)
                prefix              = match.group(1)
                parameterDictHash   = match.group(2)
                environmentDictHash = match.group(3)
                nodes               = match.group(4)
                tasks               = match.group(5)
                cores               = match.group(6)
                run                 = match.group(7)
                
                environmentDict,parameterDict,measurements = parseLikwidMetrics(resultsFolderPath + "/" + fileName,cores=="1")

                # TODO(Dominic): workaround. parameters 
                if len(environmentDict) is 0:
                   environmentDict,parameterDict,adapters = parseResultFile(resultsFolderPath + "/" + fileName.replace(".likwid",""))

                if len(measurements):
                    # write header
                    if firstFile:
                        header = []
                        header += sorted(environmentDict)
                        for parameter in knownParameters:
                            header.append(parameter)
                        for parameter in sorted(parameterDict):
                            if parameter not in knownParameters:
                                header.append(parameter)
                        header.append("nodes")
                        header.append("tasks")
                        header.append("cores")
                        header.append("run")
                        for key in sorted(measurements):
                            for subkey in measurements[key]:
                                header.append(key+" ("+subkey+")")
                        csvwriter.writerow(header)
                        firstFile=False
                    print(resultsFolderPath+"/"+fileName)
                    
                    # write row
                    row=[]
                    for key in sorted(environmentDict):
                        row.append(environmentDict[key])
                    for key in knownParameters:
                        row.append(parameterDict[key])
                    for key in sorted(parameterDict):
                        if key not in knownParameters:
                            row.append(parameterDict[key])
                    row.append(nodes)
                    row.append(tasks)
                    row.append(cores)
                    row.append(run)
                    for key in sorted(measurements):
                        for subkey in measurements[key]:
                            row.append(measurements[key][subkey])
                    csvwriter.writerow(row)

        success = not firstFile
        if success:
          # reopen the table and sort it
          tableFile   = open(tablePath, 'r')
          header      = next(tableFile)
          header      = header.strip()
          reader      = csv.reader(tableFile,delimiter=',')
          
          sortedData = sorted(reader,key=getLikwidMetricsSortingKey)
          tableFile.close()
          
          with open(tablePath, 'w') as sortedTableFile:
              writer = csv.writer(sortedTableFile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
              writer.writerow(header.split(','))
              writer.writerows(sortedData)
          print("created table:")
          print(tablePath) 

    except IOError as err:
        print ("ERROR: could not write file "+tablePath+". Error message: "<<str(err))

if __name__ == "__main__":
    import sys,os
    import configparser
    import subprocess
    import itertools
    import hashlib
    import json
    import codecs
    
    import re
    import csv
    
    subprograms = ["build","buildMissing","scripts","submit","cancel","parseAdapters","parseMetrics","cleanBuild", "cleanScripts","cleanResults","cleanAll"]
    
    knownParameters   = ["architecture", "optimisation", "dimension", "order" ]
    
    metrics =  [
                ["  MFLOP/s",                   "Sum"],  # Two whitespaces are required to not find the AVX MFLOP/s by accident
                ["AVX MFLOP/s",                 "Sum"],
                ["Memory bandwidth [MBytes/s]", "Sum"],
                ["Memory data volume [GBytes]", "Sum"],
                ["L3 bandwidth [MBytes/s]",     "Sum"], 
                ["L3 data volume [GBytes]",     "Sum"],
                ["L3 request rate",             "Avg"],
                ["L3 miss rate",                "Avg"],
                ["L2 request rate",             "Avg"],
                ["L2 miss rate",                "Avg"],
                ["Branch misprediction rate",   "Avg"]
               ]
             
    counters = [
                ["FP_ARITH_INST_RETIRED_128B_PACKED_DOUBLE", "Sum"],
                ["FP_ARITH_INST_RETIRED_SCALAR_DOUBLE",      "Sum"],
                ["FP_ARITH_INST_RETIRED_256B_PACKED_DOUBLE", "Sum"]
               ]
    
    if haveToPrintHelpMessage(sys.argv):
        info = \
"""sweep.py:

run:

./sweep.py myoptions.ini <subprogram>

available subprograms:

* build         - build all executables
* buildMissing  - build only missing executables
* scripts       - submit the generated jobs
* cancel        - cancel the submitted jobs
* parseAdapters - read the job output and parse adapter times
* parseMetrics  - read the job output and parse likwid metrics
* cleanAll      - remove the whole sweep benchmark suite
* cleanBuild    - remove the build subfolder
* cleanScripts  - remove the scripts subfolder
* cleanResults  - remove the results subfolder
* cleanHistory  - clean the submission history

typical workflow:

./sweep.py myoptions.ini build
./sweep.py myoptions.ini scripts
./sweep.py myoptions.ini submit

(after jobs have finished)

./sweep.py myoptions.ini parseAdapters

"""
        print(info)
        sys.exit()
    
    configFile = parseArgument(sys.argv,1)
    subprogram = parseArgument(sys.argv,2)
    
    config = configparser.ConfigParser()
    config.optionxform=str
    config.read(configFile)
    
    general          = config["general"]
    jobs             = config["jobs"]
    environmentSpace = parseEnvironment(config)
    parameterSpace   = parseParameters(config)
     
    exahypeRoot      = general["exahype_root"]
    outputPath       = general["output_path"]
    projectPath      = general["project_path"]
    projectName      = general["project_name"]
   
    buildFolder       = "build" 
    buildFolderPath   = exahypeRoot+"/"+outputPath+"/build"
    scriptsFolderPath = exahypeRoot+"/"+outputPath+"/scripts"
    resultsFolderPath = exahypeRoot+"/"+outputPath+"/results"
    historyFolderPath = exahypeRoot+"/"+outputPath+"/history"
    
    nodeCounts = [x.strip() for x in jobs["nodes"].split(",")]
    taskCounts = [x.strip() for x in jobs["tasks"].split(",")]
    coreCounts = [x.strip() for x in jobs["cores"].split(",")]
    runNumbers = [x.strip() for x in jobs["run"].split(",")]
    
    verifySweepAgreesWithHistoricExperiments()
    
    # select subprogram
    if subprogram == "clean":
        clean(general)
    elif subprogram == "cleanBuild":
        clean("build")
    elif subprogram == "cleanScripts":
        clean("scripts")
    elif subprogram == "cleanResults":
      clean("results")
    elif subprogram == "cleanHistory":
      clean("history")
    elif subprogram == "build":
        build()
    elif subprogram == "buildMissing":
        build(True)
    elif subprogram == "scripts":
        generateScripts()
    elif subprogram == "submit":
        submitJobs()
    elif subprogram == "cancel":
        cancelJobs()
    elif subprogram == "parseAdapters":
        parseAdapterTimes()
    elif subprogram == "parseMetrics":
        parseMetrics()
