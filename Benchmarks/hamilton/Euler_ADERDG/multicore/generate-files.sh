#!/bin/bash
#
# Perform multicore speedup tests on Hamilton.
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
> multicore/submit-all-jobs.sh
chmod u+x multicore/submit-all-jobs.sh

project=Euler_ADERDG

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
#hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
hMax=(0.0404 0.012784810126582278 0.004190871369294606 0.0013892709766162312 0.0004622425629290618) # 1/(3^l-2) times 1.01
io=no-output # or output

kernels=gen # this is just an identifier; actual kernels must be chosen before building the executables

# Derived options
i=0
mesh=regular-$i
h=${hMax[i]}

for fuseAlgorithmicSteps in 'on' 'off'
do
  for spawnPredictorAsBackgroundThread in 'on' 'off'
  do
    prefix=$project-$io-$kernels
    if [ "$fuseAlgorithmicSteps" == 'on' ]; then
      prefix+=-fused
    else
      prefix+=-nonfused
    fi
    if [ "$spawnPredictorAsBackgroundThread" == 'on' ]; then
      prefix+=+bon
    else
      prefix+=+boff
    fi
    prefix+=-$mesh

    for sharedMem in 'TBB' 'OMP'
    do       
      for order in 3 5 7
      do 
        # Create script
        script=multicore/hamilton.slurm-script
        newScript=multicore/hamilton-$prefix-$sharedMem-p$order-n1-t1.slurm-script
        cp $script $newScript
        echo "sbatch $newScript" >> multicore/submit-all-jobs.sh
       
        sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript

        sed -i -r 's,sharedMem=(TBB|OMP),sharedMem='$sharedMem',g' $newScript 
        
        sed -i 's,-p3,-p'$order',g' $newScript

        sed -i 's,script=multicore/hamilton.slurm-script,script='$newScript',g' $newScript 
  

        # Create spec files
        #for coresPerTask in 1 12 24
        for coresPerTask in 1 2 4 6 8 12 16 18 24
        #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 32 # ham6
        do
          spec=multicore/Euler_ADERDG-$io.exahype
          filename=multicore/$prefix-p$order-t1-c$coresPerTask 
          newSpec=$filename'.exahype'

          cp $spec $newSpec

          sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
          sed -i -r 's,cores(\s*)=(\s*)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
          if [ "$spawnPredictorAsBackgroundThread" == 'on' ]; then
            sed -i -r 's,configure(\s*)=(\s*)(\{.*\}),configure\1=\2{background-tasks:'$coresPerTask'},' $newSpec
          else 
            sed -i -r 's,configure(\s*)=(\s*)(\{.*\}),configure\1=\2{background-tasks:-2},' $newSpec
          fi
   
          sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
          sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
          sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
          sed -i -r 's,spawn-predictor-as-background-thread(\s*)=(\s*)(on|off),spawn-predictor-as-background-thread\1=\2'$spawnPredictorAsBackgroundThread',' $newSpec
 
          sed -i -r 's,order(\s*)const(\s*)=(\s*)([0-9]+),order\1const\2=\3'$order',' $newSpec
          sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',' $newSpec
        done
      done
    done
  done
done
