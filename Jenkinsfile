#!/usr/bin/env groovy
import hudson.FilePath;
import jenkins.model.Jenkins

def slurmBatch(code) {
sh '''#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Kill all called processes when script finishes!
trap "kill 0" SIGINT

mkdir -p ${SCRATCH}/jenkins/exahype/tmp
tmpfile=$(mktemp $SCRATCH/jenkins/exahype/tmp/slurm-XXXXX.sh)

workspace="$(pwd)"

# Write batch job header.
cat > "$tmpfile" <<EOF
#!/usr/bin/env bash
#SBATCH -o ${workspace}/job.%j.out
#SBATCH -D ${workspace}/
#SBATCH --get-user-env
#SBATCH -J exahype
#SBATCH --clusters=tum
#SBATCH --partition=acc
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
source /etc/profile.d/modules.sh
cd "${workspace}"
export OMP_NUM_THREADS=32
EOF
# Now write actual code to file.
# Quoting the EOL avoids string interpolation!
cat >> $tmpfile <<'EOF'
''' + "$code" +
'''
EOF

batch_res="$(sbatch ${tmpfile})"
echo $batch_res
jobid="$(echo ${batch_res} | grep -E -o '[0-9]*' | head -n1)"
echo "Jobid=${jobid}"

# Follow log file from batch job - gets killed after job is finished!
log_file="${workspace}/job.${jobid}.out"
touch "$log_file"
tail -f "$log_file" &

# Wait until batch job is finished...
while 
    sleep 60s # Busy waiting, sue me.
    scontrol_err="$(scontrol show job ${jobid} 2>&1 > /dev/null)"
    if grep -q Invalid <<< "${scontrol_err}"; then
	echo "Job-Id invalid/Job cancelled. Aborting."
	status="CANCELLED"
    else
	status="$(scontrol show job ${jobid} | grep -o -E "JobState=(PENDING|RUNNING|COMPLETED|FAILED|CANCELLED)" | cut -d '=' -f2)"
    fi
    [[ "$status" != "CANCELLED" && "$status" != "COMPLETED" && "$status" != "FAILED" ]]
do
continue
done

echo "status=$status"

# Kill tail.
trap 'kill $(jobs -p)' EXIT

# Delete tmp file
rm -rf "$tmpfile"

# Return correct exit code for jenkins.
if [[ "$status" == "COMPLETED" ]]; then
exit 0
else
exit 1
fi
'''
}


def build(path, name, workspace) {
    node ('mac-intel') {
	ws("${env.JOB_NAME}-${name}-build") {
slurmBatch '''
set -x
IFS=$'\n\t'
pwd
module load git subversion java/1.8 scons >/dev/null 2>&1
module unload pythonLib intel >/dev/null 2>&1
module load intel/17.0 >/dev/null 2>&1
module unload mpi.intel >/dev/null 2>&1
module load mpi.intel/2017 >/dev/null 2>&1
module unload gcc >/dev/null 2>&1
module load gcc/7 >/dev/null 2>&1
module unload tbb >/dev/null 2>&1
module load tbb/2017 cmake binutils >/dev/null 2>&1
module switch python/3.5_intel >/dev/null 2>&1

set -euo pipefail
''' + """
cp -r ${workspace}/* .
path=${path}
""" + '''
ls

# Fix linker flags
export PROJECT_LFLAGS="-lrt" 

# Build toolkit
Peano/checkout-update-peano.sh
Toolkit/build.sh

# Build example
echo "Building $path"
dir="$(readlink -f $(dirname ${path}))"

java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ${path}
cd $dir
make -j 32
'''
	    deleteDir()
	}
    }

}

def assembleBuildMatrix(combinations, workspace) {
    def buildMatrix = combinations.collectEntries {
	["${it[0]}" : {
	    -> build(it[1], it[0], "${workspace}")
	}]
    }
    return buildMatrix
}


pipeline {
    agent { label 'mac-intel' }
    // triggers {
    // 	cron('H H * * *')
    // }
    stages {
	stage ('Checkout') {
	    steps {
		sh 'echo Running as $(id -un)@$(hostname -f)'
		checkout scm
		sh 'ls'
		sh 'pwd'
		// script {
		// config = load "scripts/Jenkins/Config.groovy"
		// ...
		// }
	    }
	}

	stage ('Build') {
	    steps {
		sh 'pwd'
		sh 'ls'
		sh 'ls ApplicationExamples || true'
		script {
		    def files = findFiles(glob: 'ApplicationExamples/TestCases/*/*.exahype')
		    files.each { d -> echo d.getPath() }
		    testCases = []
		    for (file in files) {
			testCases << [file.getName(), file.getPath()]
		    }
		    def buildMatrix = assembleBuildMatrix(testCases, "${workspace}")
		    parallel buildMatrix
		}
	    }
	}
    }

	post {
	    always {
		cleanWs()
	    }
	    // success {
	    //     updateGitlabCommitStatus(name: 'jenkins', state: 'success')
	    // }
	    // failure {
	    //     updateGitlabCommitStatus(name: 'jenkins', state: 'failed')
	    // }
	}
	//options {
	//gitLabConnection('GitLab-Api-Here!')
	//}
}
