#!/bin/bash
#
# Perform multicore speedup tests on Hamilton.
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
project=Euler_FV

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
io=no-output # or output

kernels=gengodunov # this is just an identifier; actual kernels must be chosen before building the executables # gengodunov or genmusclhancock

i=0
mesh=regular-$i
h=${hMax[i]}

# Derived options
for fuseAlgorithmicSteps in "on" "off"
do

prefix=$project-$io-$kernels
if [ "$fuseAlgorithmicSteps" == "on" ]; then
  prefix+="-fused"
else
  prefix+="-nonfused"
fi
prefix+="-$mesh"

for patchSize in 7 11 15 19 # corresponds to orders=3 5 7 9
do
  # SIMULATION END TIME
  T=(  0.00333333333334 0.00111111111112 0.00037037037038 0.000123456790124 0.0000411522633746  )        # N=7
  if (( patchSize == 11 )); then
    T=( 0.00212121212122 0.000707070707072 0.0002356902357 0.0000785634118968 0.0000261878039657 )       # N=11; 11/7*T_3 ceiled with sig. 1e-6
  fi
  if (( patchSize == 15 )); then
    T=( 0.0015555555555587 0.0005185185185227 0.0001728395061773 0.0000576131687245 0.0000192043895748 ) # N=15
  fi
  if (( patchSize == 17 )); then
    T=(0.0012280701754411 0.0004093567251495 0.0001364522417189 0.0000454840805720 0.0000151613601906 )  # N=19
  fi
  t=${T[i]}
  
  for coresPerTask in 1 14 28
  do
    # Create script
    script=multicore/supermuc.load-leveler
    newScript=multicore/supermuc-$prefix-p$patchSize-n1-t1.load-leveler
    cp $script $newScript
    
    sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript
    
    sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
    sed -i 's,p3,p'$patchSize',g' $newScript
    sed -i 's,script=multicore/supermuc.load-leveler,script='$newScript',g' $newScript
  
    # Create spec file
    spec=multicore/$project-$io.exahype
    filename=multicore/$prefix-p$patchSize-t1-c$coresPerTask
    newSpec=$filename'.exahype'

    cp $spec $newSpec

    sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
    sed -i -r 's,ranks-per-node:([0-9]+),ranks-per-node:1,' $newSpec 
    sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
   
    sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
    sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
    sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
  
    sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$patchSize',' $newSpec
    sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',' $newSpec
  done
done
done
