#!/usr/bin/env groovy
util = load "Miscellaneous/Jenkins/Util.groovy"

def build(config, workspace) {
    node ('Linux-Cluster') {
	ws("${env.JOB_NAME}-${config.name}-build") {
	    deleteDir()
	    def directory=sh(returnStdout: true,
			     script: "dirname ${config.exahypeFile}"
	    ).trim()
	    echo "Directory is $directory"
util.slurmBatch '''
set -x
IFS=$'\n\t'
pwd
''' + util.getModuleCode() + """
set -euo pipefail
cp -rf ${workspace}/. .
path=${config.exahypeFile}
""" + '''
ls

# Fix linker flags
export COMPILER_LFLAGS='-pthread'
export PROJECT_LFLAGS="-lrt" 


# Build example
echo "Building $path"
dir="$(readlink -f $(dirname ${path}))"

./Toolkit2/toolkit.sh -s -j -d ${path}

cd $dir
make -j 32
'''
	    // Stash files for later reuse
	    // This is needed because the run step could run on another node!
	    // We need to exclude any example spec file from Toolit2 as these are erroneously considered in the run process
	    stash includes: "${directory}/**", excludes: "${directory}/Toolkit2/examples/*", name: "exahype-${config.name}" 
	    deleteDir()}
	    
    }

}

def assembleBuildMatrix(combinations, workspace) {
    def buildMatrix = combinations.collectEntries {
	["${it.name}" : {
	    -> build(it, "${workspace}")
	}]
    }
    return buildMatrix
}


return this
