#!/bin/bash
#
# Creates output folders in the current directory
# according t the paths given in the specfiles found in
# the arguments.
#
# Use this in order to circumvent that most ExaHyPE plotters
# don't create output directories.
#

set -e
for file in $@; do
	echo "Inspecting ExaHyPE specfile $file ..."
	for output in $(grep -A5 "plot" "$file" | grep "output" | cut -d'=' -f2); do
		d=$(dirname "$output")
		echo "mkdir $d"
		mkdir -p $(dirname "$output");
	done;
done

