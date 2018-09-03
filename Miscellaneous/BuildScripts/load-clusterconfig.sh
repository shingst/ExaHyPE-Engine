#!/bin/bash
#
# Just a neat interface to the cluster-configs directory.
# Usage is like without arguments if you have a hostname file.
# In any case, you should source this script. 
#
# The advantage of this script instead of sourcing directly the cluster
# configuration file is that it can be invoked from any directory.
#
# Usage examples:
#
#  source ..../load-clusterconfig.sh  supermuc-mpi
#  source ..../load-clusterconfig.sh  loewe
#  source ..../load-clusterconfig.sh  
#
# Without argument, will look up ~/.hostname
#
#

## remove when sourced
## if ! [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
##	>&2 echo "This script has to be sourced."
##fi

SCRIPTDIR=$(dirname "${BASH_SOURCE[0]}")

HOST_INFO_FILE="$HOME/.hostname"
CLUSTERCONFIG_DIR="$SCRIPTDIR/../ClusterConfigs"

if ! [ -z "$1" ]; then
	CLUSTERNAME="$1"
elif [[ -e "$HOST_INFO_FILE" ]]; then
        CLUSTERNAME="$(< $HOST_INFO_FILE)"
else
	>&2 echo "Usage: source $0 <NameOfCluster>"
	>&2 echo "Alternatively, put a name to $HOST_INFO_FILE"
	>&2 echo "Without sourcing, you can inspect a cluster configuration."
	>&2 echo "Available cluster configurations in $CLUSTERCONFIG_DIR :"
	>&2 echo
	for x in $(ls $CLUSTERCONFIG_DIR); do
		>&2 echo ${x%.*}
	done
	exit -1 ## remove when sourced
fi

CLUSTERCONFIG="${CLUSTERNAME}.cfg"
if [[ -e $CLUSTERCONFIG_DIR/$CLUSTERCONFIG ]]; then
	echo "# Cluster configuration for $CLUSTERNAME"
	echo "# Load this file by calling:"
	echo "#     source <(exa config $CLUSTERNAME)"
	echo "# or  source $CLUSTERCONFIG_DIR/$CLUSTERCONFIG"
	echo
	echo "echo 'Loading Cluster configuration for $CLUSTERNAME'"
	echo
	
	echo '# Clusterconfigs are supposed to be sourced from the ClusterConfig directory'
	echo 'OLDPWD=$PWD'
	echo "cd $CLUSTERCONFIG_DIR"
	
	##source $CLUSTERCONFIG
	cat $CLUSTERCONFIG_DIR/$CLUSTERCONFIG
	
	echo
	echo 'has() { type $@ &>/dev/null; } # a way to check if command is available'
	echo 'has module && module list'
	
	echo 'cd $OLDPWD'
else
	echo "## Cluster $CLUSTERNAME detected, but no ClusterConfig present!"
	echo "echo 'Cluster $CLUSTERNAME detected, but no ClusterConfig present!'"
	# echo "Create one at $CLUSTERCONFIG if you like."
fi


