#!/bin/bash
# Mandatory parameters are:
# time, nodes, tasks,
# job_name, output_file, error_file, 
# job_file, spec_file, app, 
# environment, parameters
# 
# Optional parameters are:
# ranks, cores, mail

#SBATCH -J ExaHyPE
#SBATCH -o %j_result.out
#SBATCH -t 01:59:00
#SBATCH --get-user-env
#SBATCH --exclusive
#SBATCH --cluster=mpp2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=jean-matthieu.gallard@tum.de
#SBATCH --mail-type=all

source /etc/profile.d/modules.sh

module load tbb/2017

# pipe some information into output file
echo "Timestamp (YYYY/MM/dd:hh:mm:ss): `date +%Y/%m/%d:%H:%M:%S`"
echo ""
module list
echo ""
printenv
echo ""

# multiple ranks
./ExaHyPE-Euler _spec_Euler_Limiter.exahype

echo ""
echo "Timestamp (YYYY/MM/dd:hh:mm:ss): `date +%Y/%m/%d:%H:%M:%S`"
echo ""
