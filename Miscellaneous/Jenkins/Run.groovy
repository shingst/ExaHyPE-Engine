#!/usr/bin/env groovy
util = load "Miscellaneous/Jenkins/Util.groovy"

def run(config, workspace) {
    node ('Linux-Cluster') {
	ws("${env.JOB_NAME}-${config.name}-run") {
	    unstash "exahype-${config.name}"
	    util.slurmBatch '''
set -x
ulimit -s unlimited
IFS=$'\n\t'
pwd
''' + util.getModuleCode() + '''
set -euo pipefail

ls

#file="$(find . -name '*\\.exahype2')"
#echo "File is called ${file}."
#project_name="$(grep -o "exahype-project .*" "${file}" | sed "s/exahype-project\\s*//g")"
#echo "Project is called ${project_name}"
#filename=$(basename -- "${file}")
filename=${config.name}2
#directory="$(dirname "${file}")"
directory=${config.path}
cd "${directory}"

ls
cat "${filename}"
eval "./ExaHyPE-${project_name} "${filename}""
'''
	    deleteDir()
	}
    }
}

def assembleRunMatrix(combinations, workspace) {
    def runMatrix = combinations.collectEntries {
	["${it.name}" : {
	    -> run(it, "${workspace}")
	}]
    }
    return runMatrix
}

return this
