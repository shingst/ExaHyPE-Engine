inputFiles=("./Output_baseline_ucx_patched/results/Elastic-timestep-times.csv" "./Output_tasksharing_progression_ucx_patched/results/Elastic-timestep-times.csv" "Output_sharing_server_bf_patched_ucx/results/Elastic-timestep-times.csv" "Output_tasksharing_no_progression_ucx_patched/results/Elastic-timestep-times.csv")
outputFiles=("./plots/data_2_r_2_nodes_base"  "./plots/data_2_r_2_nodes_progression" "./plots/data_2_r_2_nodes_bluefield" "./plots/data_2_r_2_nodes_no_progression" ) 

orders=( 7 8 9 )

let i=0

for table in "${inputFiles[@]}"
do
 
    echo $table
    if [ -e $table ]; then
      for o in "${orders[@]}"
      do 
        (./teaMPITaskSharing/tableslicer.py $table \
          --cols order cores consumerTasks ranks nodes replication normalised_realtime_min\
          --filter ranks=2 order=$o\
          --output "${outputFiles[i]}_$o.csv")
      done
    fi
  let i=i+1
done


