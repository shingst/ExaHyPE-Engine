#!/bin/bash

unset HOST
unset HOSTNAME

printenv

echo $TMPI_FILE

HOST=$(hostname) 
#| cut -d"." -f1
HOST=$(echo $HOST | cut -d"." -f1)
echo $HOST

echo "finished"
NODES=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODES=(${NODES})

echo $HOSTNAME

APPLICATION=$1
APP_PARAM=$2


${APPLICATION} $2 &
sleep 1
pids=($(pgrep ExaHyPE))

if [ "$HOST" = "${NODES[0]}" ]; then
  date=$(date)
  echo $HOST", "$date" , pausing ${pids[0]}"
  kill -STOP ${pids[0]}
  sleep 50
  date=$(date)
  echo $HOST", "$date" , resuming ${pids[0]}"
  kill -CONT ${pids[0]}
fi  

while( kill -0 ${pids[0]}); do
   sleep 1
done
