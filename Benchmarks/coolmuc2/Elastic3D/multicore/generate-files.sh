#!/bin/bash
#
# Perform multicore speedup tests on coolmuc2.
#
# coolmuc2 uses SLURM. SLURM supports array jobs.
#

project=Elastic3D

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
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

for order in 3 5 7 9
do
  # SIMULATION END TIME
  T=(  0.00333333333334 0.00111111111112 0.00037037037038 0.000123456790124 0.0000411522633746  )            # p=3
  if (( order == 5 )); then
    T=( 0.00212121212122 0.000707070707072 0.0002356902357 0.0000785634118968 0.0000261878039657 )  # p=5; (2*3+1)/(2*order+1)*T_3 ceiled with sig. 1e-6
  fi
  if (( order == 7 )); then
    T=( 0.0015555555555587 0.0005185185185227 0.0001728395061773 0.0000576131687245 0.0000192043895748 ) # p=7
  fi
  if (( order == 9 )); then
    T=(0.0012280701754411 0.0004093567251495 0.0001364522417189 0.0000454840805720 0.0000151613601906 ) # p=9
  fi
  t=${T[i]}
  
  # Create script
  script=multicore/coolmuc2.slurm-script
  newScript=multicore/coolmuc2-$prefix-p$order-n1-t1.slurm-script
  cp $script $newScript
 
  sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript

  sed -i 's,p3,p'$order',g' $newScript

  sed -i 's,script=multicore/coolmuc2.slurm-script,script='$newScript',g' $newScript 
  
  # Create spec files
  #for coresPerTask in 1 12 24
  for coresPerTask in 1 2 4 7 8 14 16 28
  do
    spec=multicore/Elastic3D-$io.exahype
    filename=multicore/$prefix-p$order-t1-c$coresPerTask 
    newSpec=$filename'.exahype'

    cp $spec $newSpec

    sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
    sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
    sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
   
    sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
    sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
    sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
  
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',' $newSpec
    sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',' $newSpec
  done
done
done
