#!/bin/bash

output_folder=$1
filename=$2
log_max_offload=$3

#ranks=2

cd ${output_folder}
echo "1,0,0" > ${filename}_0.txt

for i in $(seq 0 ${log_max_offload}) 
do 
j=$((2 **i))
ip=$((i+1))
echo "1,0,${j}" > ${filename}_${ip}.txt
done



