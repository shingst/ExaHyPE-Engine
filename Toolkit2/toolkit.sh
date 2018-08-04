#!/bin/bash
#
# This is a neat frontend to the ExaHyPE toolkit NG.
# It checks whether all Python is available and then calls it.
# 

Toolkit2="$(dirname $0)"
has() { type $@ &>/dev/null; } # a way to check if command is available

# find suitable python on the path
if has python3; then PYTHON3="python3";
elif python --version | grep -qi "python 3"; then PYTHON3="python"
else echo "$0: Python3 required for running the ExaHyPE toolkit" >&2; exit -1; fi

# check that all required modules are there
for module in jinja2; do
	if ! $PYTHON3 -c "import $module" 2>&1 >/dev/null; then
		echo "$0: Required python3 module '$module', not available." >&2
		exit -1
	fi
done

exec $PYTHON3 $Toolkit2/exahype/toolkit/frontend.py $@
