#!/usr/bin/env groovy
util = load "Miscellaneous/Jenkins/Util.groovy"

def run(config, workspace) {
    node ('Linux-Cluster') {
	ws(util.adjust_ws("${env.JOB_NAME}-${config.name}-run")) {
	    unstash "exahype-${config.name}"

	    def directory=sh(returnStdout: true,
			     script: "dirname ${config.exahypeFile}").trim()
	    util.slurmBatch('''
set -x
ulimit -s unlimited
IFS=$'\n\t'
pwd
''' + util.getModuleCode() + '''
set -euo pipefail

ls

file="$(find . -name '*\\.exahype2')"
echo "File is called ${file}."
#project_name="$(grep -o "exahype-project .*" "${file}" | sed "s/exahype-project\\s*//g")" //Toolkit1
project_name="$(cat ${file} | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["project_name"])')"

echo "Project is called ${project_name}"
filename=$(basename -- "${file}")
directory="$(dirname "${file}")"
cd "${directory}"

ls
cat "${filename}"
eval "mpiexec -n ${JENKINS_MPI_RANKS} ./ExaHyPE-${project_name} "${filename}""
''', directory)
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
