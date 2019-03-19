#!/bin/bash

if [ "$#" -eq 0 ]; then
	>&2 echo "Please provide path to project folder(s)"
	exit 1
fi


if [ "$#" -eq 1 ]; then
        >&2 echo "Please provide name of release build"
        exit 1
fi


folders=(ExaHyPE Toolkit CodeGenerator Peano $1)
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


exaVersion="ExaHyPE-"$2"-$(date +%Y-%m-%d)-$(git log --format="%h" -n 1)"
tarName="$exaVersion.tar.gz"

echo $tarName 

tar --exclude-vcs --exclude=*.o -czvf  $tarName ${folders[*]} ${files[*]}

