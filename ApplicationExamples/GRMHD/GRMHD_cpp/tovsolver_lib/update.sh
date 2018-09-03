#!/bin/bash
#
# This is a small script to copy an ExternalLibrary source code to the local
# application. We track changes to these codes at the ExternalLibrary directory
# in the Astrophysics repository but would like to provide these codes without
# any linking at the Main Engine repository.
# SvenK, 2018-08-29

ExternalLibraries="../../../../ExternalLibraries"

[[ -e $ExternalLibraries ]] || { echo "ExternalLibraries not found"; exit -1; }

TOVSolver="$ExternalLibraries/TOVSolver"

for remotef in $TOVSolver/src/*.{cpp,h}; do
   localf=$(basename $remotef)
   if [[ -e $localf ]] && ! diff $localf $remotef ; then
       echo "$localf has local changes, compared to $remotef , please merge back!"
       exit -1
   else
       cp -v $remotef $localf
   fi
done

echo "Done"
