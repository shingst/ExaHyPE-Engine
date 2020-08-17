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

sleep 100

while(true); do

   if kill -0 ${pids[0]} ; then
     if [ "$HOST" = "${NODES[0]}" ]; then
      #for i in ${pids[@]};
      #do
       date=$(date)
       echo $date" , sending SIGUSR to ${pids[0]}"
       kill -USR1 ${pids[0]}
      #done
     fi
   else
    exit 0
   fi
 
   sleep 10 
done
