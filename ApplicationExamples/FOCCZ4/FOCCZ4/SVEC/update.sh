#!/bin/bash
#
# This is a small script to copy the SVEC source code to the local application.
# We actually manage the SVEC code at https://bitbucket.org/svek/svec
# and want to keep the code there. ExaHyPE lacks a proper way of managing external
# repos...
# SvenK, 2018-07-17

SVEC="../../../ExternalLibraries/SVEC"

[[ -e $SVEC ]] || { echo "Cannot find $SVEC"; exit -1; }

svec="$SVEC/svec"

[[ -e $svec ]] || $SVEC/download-SVEC.sh || { echo "Cannot download SVEC"; exit -1; }

for remotef in $svec/src/PDE/tensish.cpph; do
   localf=$(basename $remotef)
   if [[ -e $localf ]] && ! diff $localf $remotef ; then
       echo "$localf has local changes, compared to $remotef , please merge back!"
       exit -1
   else
       cp -v $remotef $localf
   fi
done

echo "Done"
