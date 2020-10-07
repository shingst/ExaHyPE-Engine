#!/bin/bash

NNODES=1642
RANKS=39370

TEAMS=2
RANKS_PER_TEAM=$( echo "$RANKS/$TEAMS"  | bc )
echo "Ranks per team:"$RANKS_PER_TEAM
RANKS_PER_NODE=$( echo "($RANKS-$TEAMS+$NNODES-1)/($NNODES-1)"  | bc )
echo "Ranks per node:"$RANKS_PER_NODE
NODES_PER_TEAM=$( echo "(${RANKS_PER_TEAM}-1)/${RANKS_PER_NODE}" | bc ) 
echo "Nodes per team:"$NODES_PER_TEAM

REM_RANKS=$( echo "($RANKS-$TEAMS)-(${RANKS_PER_NODE}*${NNODES-1})" | bc )
echo "Remaining ranks:"${REM_RANKS}

NODE_LIST=()

for i in $(seq 1 ${NNODES})
do
 NODE_LIST+=("node_"${i})
done

echo ${NODE_LIST[0]}
echo ${NODE_LIST[1]}

NODE_LIST_COMP=(${NODE_LIST[@]:1:$NNODES})

MPI_STRING=""
APP_STRING=${1}
echo ${APP_STRING}

let NEXT_NODE_IDX=0
let CURR_RANKS_ON_NODE=0
let TOTAL_RANKS=0

for t in $(seq 1 $TEAMS)
do
  MPI_STRING+=" -host ${NODE_LIST[0]} -np 1 ${APP_STRING} "
  let TOTAL_RANKS+=1
  for i in $(seq 1 ${RANKS_PER_TEAM})
  do
    #for r in $(seq 1 ${RANKS_PER_NODE})
    #do
    if [ $CURR_RANKS_ON_NODE -lt $RANKS_PER_NODE ]
    then
      MPI_STRING+=" \n "
      MPI_STRING+=" -host ${NODE_LIST_COMP[${NEXT_NODE_IDX}]} -n 1 ${APP_STRING} "  
      let TOTAL_RANKS+=1
      let CURR_RANKS_ON_NODE+=1
    else
      let NEXT_NODE_IDX=NEXT_NODE_IDX+1
      let CURR_RANKS_ON_NODE=0
      MPI_STRING+=" \n "
      MPI_STRING+=" -host ${NODE_LIST_COMP[${NEXT_NODE_IDX}]} -n 1 ${APP_STRING} "  
      let TOTAL_RANKS+=1
      let CURR_RANKS_ON_NODE+=1
    fi
  done
done

echo "Total ranks in this setup:"${TOTAL_RANKS}

echo -e ${MPI_STRING} > test.out
