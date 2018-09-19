#!/bin/bash

TOOLKIT="java -jar ./Toolkit/dist/ExaHyPE.jar --not-interactive"

# code output directories
rm -rf output
mkdir -p output/test-constants-{aderdg,fv,limiting}

examples="Toolkit/examples"

cd  ../../

for specfile in $examples/test-constants-{aderdg,fv,limiting}.exahype; do
#for specfile in $examples/test-constants-limiting.exahype; do
	outputdir=$(basename "$specfile")
	outputdir=$examples/output/${outputdir%.*}
	
	(
	$TOOLKIT $specfile && cd $outputdir && make && echo "SUCCESS with $specfile" || echo "FAILURE with $specfile"
	) 
done
