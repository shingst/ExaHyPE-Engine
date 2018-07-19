#!/bin/bash

#from: https://stackoverflow.com/questions/24112727/relative-paths-based-on-file-location-instead-of-current-working-directory#
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

if [ $# -eq 0 ]; then
    if [ ! -f ${parent_path}/../Submodules/Peano/.git ]; then
	echo "Initialize Peano submodule"
	git submodule update --init --remote
    else
       echo "Update Peano submodule"
       cd ${parent_path}/../Submodules/Peano
       git pull origin master
    fi
else
    while getopts hsw opt; do
	case $opt in
	    h) echo "-h prints this message"
	       echo "-s set Peano Submodule url to ssh"
	       echo "-v set Peano Submodule url to https";;
	    s) git config submodule.Submodules/Peano.url git@gitlab.lrz.de:gi26det/Peano.git;;
	    w) git config submodule.Submodules/Peano.url https://gitlab.lrz.de/gi26det/Peano.git;;
	esac
    done
	  
fi






