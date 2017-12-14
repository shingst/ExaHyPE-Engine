#!/bin/bash
#
# Perform multicore speedup tests on Coolmuc.
#
# Coolmuc uses SLURM. SLURM supports array jobs.
#

# PREAMBLE
project=CCZ4
order=3
skipReductionInBatchedTimeSteps=on
batchFactor=0.8
io=no-output # or output
kernels=opt # this is just an identifier; actual kernels must be chosen before building the executables
sharedMem=TBB

# MESH
i=0
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
#hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
hMax=( 0.12 )
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
  
  for nodes in 1
  do
    for tasksPerNode in 1 2 4 8 16 32 64 # knl
    #for tasksPerNode in 1 2 4 6 8 12 24 # ham7
    #for tasksPerNode in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 # ham6
    do 
      let tasks=$nodes*$tasksPerNode
      let coresPerTask=64/$tasksPerNode # knl
      #let coresPerTask=24/$tasksPerNode # ham7
      #let coresPerTask=16/$tasksPerNode # ham6
     
      # Create script
      script=single-node/coolmuc.slurm-script
      newScript=single-node/coolmuc-$prefix-p$order-n$nodes-t$tasksPerNode-c$coresPerTask-$sharedMem.slurm-script
      cp $script $newScript
     
      sed -i -r 's,--nodes(\s*)=(\s*)([0-9]*),--nodes\1=\2'$nodes',' $newScript
      sed -i -r 's,--ntasks-per-node(\s*)=(\s*)([0-9]*),--ntasks-per-node\1=\2'$tasksPerNode',' $newScript
      
      sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
      sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
      sed -i 's,p3,p'$order',g' $newScript
      
      sed -i 's,tasksPerNode=1,tasksPerNode='$tasksPerNode',' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',' $newScript
      
      sed -i 's,script='$script',script='$newScript',g' $newScript 
      
      # Create spec file
      spec=single-node/CCZ4-$io.exahype
      filename=single-node/$prefix-p$order-t$tasksPerNode-c$coresPerTask
      newSpec=$filename'.exahype'
      cp $spec $newSpec
      
      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
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
