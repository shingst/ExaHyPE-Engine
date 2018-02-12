#!/bin/bash
#
# Perform $directory speedup tests on Coolmuc.
#
# Coolmuc uses SLURM. SLURM supports array jobs.
#

project=CCZ4
directory=multicore

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
#hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
hMax=( 0.12 )
io=no-output # or output

kernels=opt # this is just an identifier; actual kernels must be chosen before building the executables

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do
i=0
mesh=regular-$i
h=${hMax[i]}

prefix=$project-$io-$kernels
if [ "$fuseAlgorithmicSteps" == "on" ]; then
  prefix+="-fused"
else
  prefix+="-nonfused"
fi
prefix+="-$mesh"

for order in 3 6 9
do
  # Create script
  script=$directory/coolmuc.slurm-script
  newScript=$directory/coolmuc-$prefix-p$order-n1-t1.slurm-script
  cp $script $newScript
 
  sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
  sed -i 's,-p3,-p'$order',g' $newScript
  sed -i 's,script=$directory/coolmuc.slurm-script,script='$newScript',g' $newScript 
  
  # Create spec files
  for coresPerTask in 1 2 4 8 16 32 64
  #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 48 # ham7
  #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 32 # ham6
  do
  for bg in "on" "off"
  do
    spec=$directory/$project-$io.exahype
    filename=$directory/$prefix-p$order-t1-c$coresPerTask-b${bg}
    newSpec=$filename'.exahype'

    cp $spec $newSpec
    
    sed -i -r 's,spawn-predictor-as-background-thread(\s*)=(\s*)(on|off),spawn-predictor-as-background-thread\1=\2'${bg}',' $newSpec
    sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
    sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
    if [ $bg == "on" ]; then
        sed -i -r 's,background-tasks(\s*):(\s*)([0-9]+),background-tasks\1:\2'$coresPerTask',' $newSpec
    else
        sed -i -r 's,background-tasks(\s*):(\s*)([0-9]+),background-tasks\1:\2'-2',' $newSpec
    fi
    sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
    sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
    sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
  
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',' $newSpec
    sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',' $newSpec
  done
  done
done
done