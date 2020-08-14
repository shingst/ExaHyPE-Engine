#!/bin/bash

unset HOST
unset HOSTNAME

printenv

echo $TMPI_FILE

echo "try 1:"
#hostname

echo "try 2:"
HOST=$(hostname) 
#| cut -d"." -f1
HOST=$(echo $HOST | cut -d"." -f1)
echo $HOST

echo "finished"
NODES=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODES=(${NODES})

echo $HOSTNAME

#printenv
#echo ${NODES[0]}
#echo ${NODES[1]}

APPLICATION=$1
APP_PARAM=$2

#echo $APPLICATION
#echo $APP_PARAM

${APPLICATION} $2 &
sleep 1
pids=($(pgrep ExaHyPE))

# offset when to start disturbing
sleep 100

export increment=0.1
export current_sleep=0.1

while(true); do

   if kill -0 ${pids[0]} ; then
     if [ "$HOST" = "${NODES[0]}" ]; then
      #for i in ${pids[@]};
      #do
       date=$(date)
       echo $HOST", "$date" , pausing ${pids[0]}"
       kill -STOP ${pids[0]}
       sleep $current_sleep
       current_sleep=$(echo $current_sleep+$increment | bc)
       date=$(date)
       echo $HOST", "$date" , resuming ${pids[0]}"
       kill -CONT ${pids[0]}
      #done
     fi
   else
    exit 0
   fi
 
   sleep 0.5
done
