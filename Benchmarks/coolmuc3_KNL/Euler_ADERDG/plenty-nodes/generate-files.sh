#!/bin/bash
#
# Perform with speedup tests on coolmuc3_KNL using several nodes.
#
# coolmuc3_KNL uses SLURM. SLURM supports array jobs.
#


# PREAMBLE
project=Euler_ADERDG
skipReductionInBatchedTimeSteps=on
batchFactor=0.8
io=no-output # or output
kernels=gen # this is just an identifier; actual kernels must be chosen before building the executables
sharedMem=None

# MESH
i=1
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
mesh=regular-$i
h=${hMax[i]}

for order in 3 5 7
do

for fuseAlgorithmicSteps in "on" "off"
do
  prefix=$project-$io-$kernels
  if [ "$fuseAlgorithmicSteps" == "on" ]; then
    prefix+="-fused"
  else
    prefix+="-nonfused"
  fi
  prefix+="-$mesh"

  for nodes in 29   # 3d
  #for nodes in 10 28 82 # 2d
  do
    for tasksPerNode in 1 28 
    do 
      let tasks=$nodes*$tasksPerNode
      let coresPerTask=64/$tasksPerNode # ham7
      #let coresPerTask=16/$tasksPerNode # ham6

      # Create script
      script=plenty-nodes/coolmuc3_KNL.slurm-script
      newScript=plenty-nodes/coolmuc3_KNL-$prefix-p$order-n$nodes-t$tasksPerNode-c$coresPerTask-$sharedMem.slurm-script
      cp $script $newScript
     
      sed -i -r 's,--nodes(\s*)=(\s*)([0-9]*),--nodes\1=\2'$nodes',' $newScript
      sed -i -r 's,--ntasks-per-node(\s*)=(\s*)([0-9]*),--ntasks-per-node\1=\2'$tasksPerNode',' $newScript
      sed -i -r 's,--cpus-per-task(\s*)=(\s*)([0-9]*),--cpus-per-task\1=\2'$coresPerTask',' $newScript
      
      sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
      sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
      sed -i 's,-p3,-p'$order',g' $newScript

      sed -i -r 's,nodes(\s*)=(\s*)([0-9]*),nodes\1=\2'$nodes',' $newScript
      sed -i 's,tasksPerNode=1,tasksPerNode='$tasksPerNode',' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',' $newScript

      sed -i 's,script='$script',script='$newScript',g' $newScript 

      # Create spec file
      spec=plenty-nodes/Euler_ADERDG-$io.exahype
      filename=plenty-nodes/$prefix-p$order-t$tasksPerNode-c$coresPerTask
      newSpec=$filename'.exahype'
      cp $spec $newSpec

      sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$tasksPerNode',g' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
     
      sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(off|on),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
     
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',g' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',g' $newSpec
     
      sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',g' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec
    done
  done
done
done
