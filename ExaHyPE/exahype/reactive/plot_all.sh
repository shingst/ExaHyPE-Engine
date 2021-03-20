#!/bin/bash

plotting_script=~/ExaHyPE-Engine-merge/ExaHyPE/exahype/stealing/stealing_statistics.py
folder=$1

cd $folder
files=$(ls *.out)

echo $files

for f in $files
do
  echo "plotting file:"${f}
  python ${plotting_script} $f
done


