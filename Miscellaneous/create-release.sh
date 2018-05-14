#!/bin/bash

if [ "$#" -eq 0 ]; then
	>&2 echo "Please provide path to project folder(s)"
	exit 1
fi



folders=(ExaHyPE Toolkit CodeGenerator Peano $@)
files=(LICENSE.txt)
for i in ${folders[@]}; do
	if [ ! -d "$i" ]; then
		>&2 echo "Cannot find folder: $i"
		exit 1
	fi
done

for i in ${files[@]}; do
	if [ ! -f "$i" ]; then
		>&2 echo "Cannot find file: $i"
		exit 1
	fi
done


#peanoVersion="Peano-$(git log --format="%h" -n 1 Peano/peano)"
exaVersion="ExaHyPE-$(git log --format="%h" -n 1)"
tarName="$exaVersion.tar.gz"

echo $tarName 

tar --exclude-vcs -czvf  $tarName ${folders[*]} ${files[*]}