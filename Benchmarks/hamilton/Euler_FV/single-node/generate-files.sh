#!/bin/bash
#
# Perform with speedup tests on Hamilton using several nodes.
#
# Hamilton uses SLURM. SLURM supports array jobs.
#
# System specification(s):
#
# Hamilton 6 (x122 nodes):
#    2 x Intel Xeon E5-2650 v2 (Ivy Bridge) 8 cores, 2.6 GHz processors (16 cores per node)
#    64 GB DDR3 memory (4 GB per core)
#    the nodes are diskless
#    1 x TrueScale 4 x QDR single-port InfiniBand interconnect
#
# Hamilton 7 (x112 nodes):
#   2 x Intel Xeon E5-2650 v4 (Broadwell) 12 cores, 2.2 GHz processors (24 cores per node)
#   64 GB TruDDR4 memory
#   the nodes are diskless
#   1 x Intel OmniPath 100 Gb InfiniBand interconnect

# PREAMBLE
project=Euler_FV
patchSize=7
skipReductionInBatchedTimeSteps=on
batchFactor=0.8
io=no-output # or output
kernels=gengodunov # this is just an identifier; actual kernels must be chosen before building the executables # gengodunov or genmusclhancock
sharedMem=None

# MESH
i=0
hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
mesh=regular-$i
h=${hMax[i]}

# SIMULATION END TIME patchSize=7 11 15 19 corresponds to orders=3 5 7 9
T=( 0.01 0.00334 0.00112 0.00038 0.00013 )            # p=3
if (( patchSize == 11 )); then
  T=( 0.006364 0.002126 0.000713 0.000242 0.000083 )  # p=5; (2*3+1)/(2*order+1)*T_3 ceiled with sig. 1e-6
fi
if (( patchSize == 15 )); then
  T=( 0.004667 0.001559 0.000523 0.000178 0.000061 )  # p=7
fi
if (( patchSize == 19 )); then
  T=( 0.003685 0.001231 0.000413 0.00014 0.000048 )   # p=9
fi
t=${T[i]}

for fuseAlgorithmicSteps in "on" "off"
do
  prefix=$project-$io-$kernels
  if [ "$fuseAlgorithmicSteps" == "on" ]; then
    prefix+="-fused"
  else
    prefix+="-nonfused"
  fi
  prefix+="-$mesh"

  for nodes in 1
  do
    for tasksPerNode in 1 2 4 8 12 24 # ham7
    #for tasksPerNode in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 # ham6
    do 
      let tasks=$nodes*$tasksPerNode
      let coresPerTask=24/$tasksPerNode # ham7
      #let coresPerTask=16/$tasksPerNode # ham6

      # Create script
      script=single-node/hamilton.slurm-script
      newScript=single-node/hamilton-$prefix-p$patchSize-n$nodes-t$tasksPerNode-c$coresPerTask-$sharedMem.slurm-script
      cp $script $newScript
     
      sed -i -r 's,--nodes(\s*)=(\s*)([0-9]*),--nodes\1=\2'$nodes',' $newScript
      sed -i -r 's,--ntasks-per-node(\s*)=(\s*)([0-9]*),--ntasks-per-node\1=\2'$tasksPerNode',' $newScript
      
      sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
      sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
      sed -i 's,p3,p'$patchSize',g' $newScript

      sed -i 's,nodes=1,nodes='$nodes',' $newScript
      sed -i 's,tasksPerNode=1,tasksPerNode='$tasksPerNode',' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',' $newScript

      sed -i 's,script=hamilton.slurm-script,script='$newScript',g' $newScript 

      # Create spec file
      spec=single-node/$project-$io.exahype
      filename=single-node/$prefix-p$patchSize-t$tasksPerNode-c$coresPerTask
      newSpec=$filename'.exahype'
      cp $spec $newSpec

      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
      sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$tasksPerNode',g' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
     
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',g' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',g' $newSpec
     
      sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$patchSize',' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec
    done
  done
done