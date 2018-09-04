#!/bin/bash
#SBATCH --partition=iboga
#SBATCH --time=10:00:00
#SBATCH --ntasks=200
#SBATCH --job-name=exahype
#SBATCH --mem-per-cpu=2600
#SBATCH --mail-type=ALL

set -e

source <(exa config iboga-gcc-mpi)
source <(exa config iboga-gcc-tbb)

cd $(exa root)/ApplicationExamples/GRMHD/GRMHD/

OUTFILE="run-$(date +%Y-%m-%dT%H-%M-%S).output"

# 5 threads each
exec mpirun.mpich -np 40 ./ExaHyPE-GRMHD -p ../GRMHD_ADERDG.exahype 2>&1 | tee $OUTFILE
