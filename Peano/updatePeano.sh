#!/bin/bash
#
# Peano is nowadays managed as a git submodule. That makes it easy
# to keep in sync. However, you also need access to the Peano
# repository (you find the URL in the file ../.gitmodules).
#
# If you don't have (read) access to the repository, you can use
# the script ./fetchSnapshot.sh  to download a public available
# tarball instead.
#


#from: https://stackoverflow.com/questions/24112727/relative-paths-based-on-file-location-instead-of-current-working-directory#
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

if [ $# -eq 0 ]; then
    if [ ! -f ${parent_path}/../Submodules/Peano/.git ]; then
	echo "Initialize Peano submodule"
	exec git submodule update --init --remote
    else
       echo "Update Peano submodule"
       cd ${parent_path}/../Submodules/Peano
       exec git pull origin master
    fi
else
    while getopts hsw opt; do
	case $opt in
	    h) echo "-h prints this message"
	       echo "-s set Peano Submodule url to ssh"
	       echo "-v set Peano Submodule url to https"
	       exit -1;;
	    s) exec git config submodule.Submodules/Peano.url git@gitlab.lrz.de:gi26det/Peano.git;;
	    w) exec git config submodule.Submodules/Peano.url https://gitlab.lrz.de/gi26det/Peano.git;;
	esac
    done	  
fi






