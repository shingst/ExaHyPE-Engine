#!/bin/bash

ranks=28
threads=2
timestep_interval=1000

plotting_tool=~/Codes/ExaHyPE-Engine/Miscellaneous/OffloadingAnalysis/plot_stp_trace.py

input_dir=$1
output_dir=$2
its=$3

echo ${input_dir}
echo ${output_dir}
echo ${its}

for i in `seq 1 ${its}`
do
 timestep=$((i*timestep_interval))
 echo ${timestep}
 python3 ${plotting_tool} ${input_dir} ${ranks} ${threads} ${timestep}
done


