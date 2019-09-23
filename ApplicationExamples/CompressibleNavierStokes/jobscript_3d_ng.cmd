#!/bin/bash
#SBATCH --account=pr83no

# Job Name and Files (also --job-name)
#SBATCH -J exa3DCloud
#SBATCH --partition=micro

#Output and error (also --output, --error):
#SBATCH -o /hppfs/scratch/0E/ga24dib3/exahype_bench/logs/%x.%j.out
#SBATCH -e /hppfs/scratch/0E/ga24dib3/exahype_bench/logs/%x.%j.err

#Initial working directory (also --chdir):
#SBATCH -D ./

#Notification and type
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=lukas.krenz@in.tum.de

# Wall clock limit:
#SBATCH --time=6:59:00
#SBATCH --no-requeue

# set number of nodes
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=12
#SBATCH --ear=off

#--constraint="scratch&work"

module load python/3.6_intel # toolkit
module load slurm_setup

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Try to fix bug with recv buffer
export PSM_MQ_REVCREQ_MAX=10000000

# run the application
mpiexec -np $SLURM_NTASKS /dss/dsshome1/0E/ga24dib3/src/ExaHyPE-Engine/ApplicationExamples/CompressibleNavierStokes/ExaHyPE-NavierStokes-3D /dss/dsshome1/0E/ga24dib3/src/ExaHyPE-Engine/ApplicationExamples/CompressibleNavierStokes/NavierStokesBubblesBench3D.exahype2
