#!/usr/bin/env groovy
def slurmBatch(code, directory) {
    sh '''#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Kill all called processes when script finishes!
trap "kill 0" SIGINT

# Set environment defaults
export JENKINS_MPI_RANKS=1
export JENKINS_MPI_NODES=1
export JENKINS_TBB_CORES=28
# Read env config file
''' + 'export FILENAME=' + directory + '/env.sh' + '''
echo "filename is ${FILENAME}"
[ -f "${FILENAME}" ] && cat "${FILENAME}"
set -o allexport # no need for export
[ -f "${FILENAME}" ] && source "${FILENAME}"
set +o allexport

echo "Using ${JENKINS_MPI_RANKS} ranks on ${JENKINS_MPI_NODES} nodes."

#REMARK: required until SCRATCH env got fixed by LRZ
export SCRATCH="/gpfs/scratch/pr63so/di57wuf"

mkdir -p ${SCRATCH}/jenkins/exahype/tmp
tmpfile=$(mktemp $SCRATCH/jenkins/exahype/tmp/slurm-XXXXX.sh)

workspace="$(pwd)"

echo $workspace
# Write batch job header.
cat > "$tmpfile" <<EOF
#!/usr/bin/env bash
#SBATCH -o ${workspace}/job.%j.out
#SBATCH -D ${workspace}/
#SBATCH --get-user-env
#SBATCH -J exahype
#SBATCH --clusters=mpp2
##SBATCH --partition=acc
#SBATCH --cpus-per-task=${JENKINS_TBB_CORES}
#SBATCH --time=02:00:00
#SBATCH --ntasks=${JENKINS_MPI_RANKS}
#SBATCH --nodes=${JENKINS_MPI_NODES}
source /etc/profile.d/modules.sh
cd "${workspace}"
export OMP_NUM_THREADS=${JENKINS_TBB_CORES}
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
    scontrol_err="$(scontrol --cluster=mpp2 show job ${jobid} 2>&1 > /dev/null)"
    if grep -q Invalid <<< "${scontrol_err}"; then
	echo "Job-Id invalid/Job cancelled. Aborting."
	status="CANCELLED"
    else
	status="$(scontrol --cluster=mpp2 show job ${jobid} | grep -o -E "JobState=(PENDING|RUNNING|COMPLETED|FAILED|CANCELLED)" | cut -d '=' -f2)"
    fi
    [[ "$status" != "CANCELLED" && "$status" != "COMPLETED" && "$status" != "FAILED" ]]
do
continue
done

echo "status=$status"

# Kill tail.
trap 'kill $(jobs -p)' EXIT

# Delete tmp file
#rm -rf "$tmpfile"

# Return correct exit code for jenkins.
if [[ "$status" == "COMPLETED" ]]; then
exit 0
else
exit 1
fi
'''
}

def getModuleCode() {
    return '''
source /etc/profile.d/modules.sh

module switch intel/18.0 >/dev/null 2>&1
module unload mpi.intel >/dev/null 2>&1
module load mpi.intel/2018 >/dev/null 2>&1
module unload gcc >/dev/null 2>&1
module load gcc/8 >/dev/null 2>&1
module unload tbb >/dev/null 2>&1
module load tbb/2018 >/dev/null 2>&1
module load cmake >/dev/null 2>&1
module switch python/3.5_intel >/dev/null 2>&1
module list
'''
}

def adjust_ws(dir){
    return dir.replaceAll('%2F','_')
}

return this
