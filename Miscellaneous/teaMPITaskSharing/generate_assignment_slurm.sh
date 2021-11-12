#!/bin/bash

NODES=${SLURM_NNODES}
RANKS=$1


if [ -z ${TEAMS} ]; then
  TEAMS=2
fi

RANKS_PER_TEAM=$( echo "$RANKS/$TEAMS"  | bc )
RANKS_PER_NODE=$( echo "($RANKS-$TEAMS+$NODES-1)/($NODES-1)"  | bc )

NODES_PER_TEAM=$( echo "(${RANKS_PER_TEAM}-1)/${RANKS_PER_NODE}" | bc ) 

NODE_LIST=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODE_LIST=(${NODE_LIST})


NODE_LIST_COMP=(${NODE_LIST[@]:1:$NODES})

MPI_STRING=""
APP_STRING=$2

let NEXT_NODE_IDX=0
let CURR_RANKS_ON_NODE=0

for t in $(seq 1 $TEAMS)
do
  MPI_STRING+=" -host ${NODE_LIST[0]} -np 1 ${APP_STRING} "
  for i in $(seq 1 ${RANKS_PER_TEAM})
  do
    if [ $CURR_RANKS_ON_NODE -lt $RANKS_PER_NODE ]
    then
      MPI_STRING+=" \n "
      MPI_STRING+=" -host ${NODE_LIST_COMP[${NEXT_NODE_IDX}]} -n 1 ${APP_STRING} "  
      let CURR_RANKS_ON_NODE+=1
    else 
      let NEXT_NODE_IDX=NEXT_NODE_IDX+1
      let CURR_RANKS_ON_NODE=0
      MPI_STRING+=" \n "
      MPI_STRING+=" -host ${NODE_LIST_COMP[${NEXT_NODE_IDX}]} -n 1 ${APP_STRING} "  
      let CURR_RANKS_ON_NODE+=1
    fi
  done
  if [ $t -lt $TEAMS ] 
  then
    MPI_STRING+=" \n "
  fi
done

echo  -e ${MPI_STRING}
