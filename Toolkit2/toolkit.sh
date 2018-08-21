#!/bin/bash
#
# This is a neat frontend to the ExaHyPE toolkit NG.
# It checks whether all Python is available and then calls it.
# 

# local var to resolve relative path correctly
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")

# find suitable python on the path
has() { type $@ &>/dev/null; } # a way to check if command is available
if has python3; then PYTHON3="python3";
elif python --version | grep -qi "python 3"; then PYTHON3="python"
else echo "$0: Python3 required for running the ExaHyPE toolkit" >&2; exit -1; fi

# Run program using "exec", which ensures the proper return value of the shell script
exec $PYTHON3 "$scriptDir"/exahype/toolkit $@
