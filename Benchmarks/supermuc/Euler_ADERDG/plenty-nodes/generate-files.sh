#!/bin/bash
#
# Perform with speedup tests on Hamilton using several nodes.
#
# SuperMUC Phase 2 uses LOADLEVELER. 
# LOADLEVELER does not support array jobs.
#
# System specification:
#
# SuperMUC Phase 2 (x3072 nodes):
#   2x Haswell Xeon Processor E5-2697 v3 (2x14 cores)
#   64 GB memory
#   the nodes are diskless
#   Infinityband FDR14 interconnects

# PREAMBLE
project=Euler_ADERDG
order=3
skipReductionInBatchedTimeSteps=on
batchFactor=0.8
io=no-output # or output
kernels=gen # this is just an identifier; actual kernels must be chosen before building the executables
sharedMem=None

# MESH
i=0
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
mesh=regular-$i
h=${hMax[i]}

# SIMULATION END TIME
T=( 0.01 0.00334 0.00112 0.00038 0.00013 )            # p=3
if (( order == 5 )); then
  T=( 0.006364 0.002126 0.000713 0.000242 0.000083 )  # p=5; (2*3+1)/(2*order+1)*T_3 ceiled with sig. 1e-6
fi
if (( order == 7 )); then
  T=( 0.004667 0.001559 0.000523 0.000178 0.000061 )  # p=7
fi
if (( order == 9 )); then
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

  for nodes in 10 28 82
  do
    for tasksPerNode in 1 2 4 7 14 28
    do 
      let tasks=$nodes*$tasksPerNode
      let coresPerTask=28/$tasksPerNode

      # Create script
      script=plenty-nodes/supermuc.load-leveler
      newScript=plenty-nodes/supermuc-$prefix-p$order-n$nodes-t$tasksPerNode-c$coresPerTask-$sharedMem.load-leveler
      cp $script $newScript
     
      queue=micro
      if (( nodes > 20 )); then
        queue=general
      fi
      sed -i -r 's,@(\s*)class(\s*)=(\s*)(\w+),@\1class\2=\3'$queue',' $newScript
      sed -i -r 's,@(\s*)node(\s*)=(\s*)([0-9]+),@\1node\2=\3'$nodes',' $newScript
      sed -i -r 's,@(\s*)tasks_per_node(\s*)=(\s*)([0-9]+),@\1tasks_per_node\2=\3'$tasksPerNode',' $newScript
      
      sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
      sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
      sed -i 's,p3,p'$order',g' $newScript

      sed -i 's,nodes=1,nodes='$nodes',' $newScript
      sed -i 's,tasksPerNode=1,tasksPerNode='$tasksPerNode',' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',' $newScript

      sed -i 's,script=plenty-nodes/supermuc.load-leveler,script='$newScript',g' $newScript

      # Create spec file
      spec=plenty-nodes/Euler_ADERDG-$io.exahype
      filename=plenty-nodes/$prefix-p$order-t$tasksPerNode-c$coresPerTask
      newSpec=$filename'.exahype'
      cp $spec $newSpec

      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
      sed -i -r 's,ranks-per-node:([0-9]+),ranks-per-node:'$tasksPerNode',g' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
     
      sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(off|on),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
      
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',g' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',g' $newSpec
     
      sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',g' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec
    done
  done
done
