#!/bin/bash
#
# Create three specfiles, one for pure ADERDG, one for pure FV, one for the LimitingADERDG
# scheme, run the toolkit on each, compile each. This test should just pass.
#
# This test is a demonstration of the variable injection capabilities of the mexa.py tool
# which goes quite well with scripting.

mexa="../mexa.py"
toolkit="java -jar ./Toolkit/dist/ExaHyPE.jar --not-interactive"

rm -rf output *.out-exahype # code output directories
examples="$(dirname $(readlink -f $0))" # absolute path to the example directory

for solver in DGSolver FVSolver LimitingSolver; do (
	output_dir="$examples/output/test-${solver}"
	specfile="$examples/test-${solver}.out-exahype"
	mkdir -p $output_dir

	$mexa --outformat exahype --infile generic-test.par \
		"Solver <= $solver" \
		"Project::Directory = '$output_dir'"\
		> $specfile

	(cd  ../../../ && $toolkit $specfile ) && \
	cd $output_dir && make -j5 && \
	echo "SUCCESS with Solver=$solver" || echo "FAILURE with Solver=$solver"
) & done

wait
echo "Resulting binaries:"
ls -lh output/*/ExaHyPE-*
