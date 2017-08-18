#!/bin/bash
#
# Perform convergence tests on Hamilton.
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
project=Euler_ADERDG

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
hMax=( 0.11112 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
T=0.03 # SIMULATION END TIME # chosen the same for all patch and mesh sizes

kernels=gen # this is just an identifier; actual kernels must be chosen before building the executables

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do 
  # prefix= project-kernels-mesh-i-
  # regex = .+-.+-(\w+)-.+-.+-
  for order in 3 5 7 9
  do
    # Create script
    prefix=$project-$kernels
    if [ "$fuseAlgorithmicSteps" == "on" ]; then
      prefix+="-fused"
    else
      prefix+="-nonfused"
    fi
    script=convergence/hamilton.slurm-script
    newScript=convergence/hamilton-$prefix-p${order}.slurm-script
    cp $script $newScript
   
    sed -i 's,prefix='$project',prefix='$prefix',g' $newScript
    sed -i 's,kernels=gen,kernels='$kernels',g' $newScript
    sed -i 's,p3,p'$order',g' $newScript
    sed -i 's,regular-0,'$mesh',g' $newScript
    sed -i 's,script=convergence/hamilton.slurm-script,script='$newScript',g' $newScript
  
    for i in 0 1 2
    do
      mesh=regular-$i
      h=${hMax[i]}
      t=$T
       
      specPrefix=$prefix-$mesh   
    
      # Create spec files
      coresPerTask=24
      spec=convergence/Euler_ADERDG.exahype
      filename=convergence/$specPrefix-p$order
      newSpec=$filename'.exahype'

      cp $spec $newSpec

      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
      sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
     
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
      sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
    
      sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.|-|\+|e|E)*),maximum-mesh-size\1=\2'$h',' $newSpec
    done
  done
done