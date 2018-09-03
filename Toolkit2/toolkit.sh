#!/bin/bash
#
# This is a neat frontend to the ExaHyPE toolkit NG.
# It checks whether all Python is available and then calls it.
# 

Toolkit2="$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")"
has() { type $@ &>/dev/null; } # a way to check if command is available
function join_by { local IFS="$1"; shift; echo "$*"; } # join bash-array with delimiter

# find suitable python on the path
has() { type $@ &>/dev/null; } # a way to check if command is available
if has python3; then PYTHON3="python3";
elif python --version | grep -qi "python 3"; then PYTHON3="python"
else echo "$0: Python3 required for running the ExaHyPE toolkit" >&2; exit -1; fi

# check that all required modules are there.
modules=(attr pyrsistent markupsafe jinja2 jsonschema)

if ! $PYTHON3 -c "import sys; sys.path.append(\"$Toolkit2\"); import $(join_by "," ${modules[@]})" 2>&1 >/dev/null; then
  echo "$0: At least one required Python3 module is not available." >&2
  if echo "$@" | grep -q -- '--interactive'; then
    ./install-dependencies.sh >&2 || { echo "$0: Installing dependencies failed."; exit -1; }
  else 
    echo "$0: Call with --interactive to let the toolkit install the dependencies locally and interactively." >&2
    exit -1
  fi
fi

#exec $PYTHON3 $Toolkit2/exahype/toolkit/frontend.py $@

# Run program using "exec", which ensures the proper return value of the shell script
exec $PYTHON3 "$Toolkit2"/exahype/toolkit $@
