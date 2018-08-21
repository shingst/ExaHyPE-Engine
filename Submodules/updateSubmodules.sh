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

# local var to resolve relative path correctly
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
currentLocation=$(pwd)

# move to the CodeGenerator directory
cd "$scriptDir"

# By default don't rebuild LIBXSMM is nothing changed
REBUILD_LIBXSMM=false

if [ $# -eq 0 ]; then
	#Peano
	if [ ! -d Peano ]; then
		mkdir Peano
	fi
	if [ ! -f Peano/.git ]; then
		echo "Initialize Peano submodule"
		git submodule update --init Peano
	else
		echo "Update Peano submodule"
		cd Peano
		git pull origin master
		cd ..
	fi
	#Jinja2
	if [ ! -d jinja ]; then
		mkdir jinja
	fi
	if [ ! -f jinja/.git ]; then
		echo "Initialize jinja submodule"
		git submodule update --init jinja
	else
		echo "Update jinja submodule"
		cd jinja
		git pull origin master
		cd ..
	fi
	#Markupsafe
	if [ ! -d markupsafe ]; then
		mkdir markupsafe
	fi
	if [ ! -f markupsafe/.git ]; then
		echo "Initialize markupsafe submodule"
		git submodule update --init markupsafe
	else
		echo "Update markupsafe submodule"
		cd markupsafe
		git pull origin master
		cd ..
	fi
	#attrs
	if [ ! -d attrs ]; then
		mkdir attrs
	fi
	if [ ! -f attrs/.git ]; then
		echo "Initialize attrs submodule"
		git submodule update --init attrs
	else
		echo "Update attrs submodule"
		cd attrs
		git pull origin master
		cd ..
	fi
	#pyrsistent
	if [ ! -d pyrsistent ]; then
		mkdir pyrsistent
	fi
	if [ ! -f pyrsistent/.git ]; then
		echo "Initialize pyrsistent submodule"
		git submodule update --init pyrsistent
	else
		echo "Update pyrsistent submodule"
		cd pyrsistent
		git pull origin master
		cd ..
	fi
	#pyrsistent
	if [ ! -d jsonschema ]; then
		mkdir jsonschema
	fi
	if [ ! -f jsonschema/.git ]; then
		echo "Initialize jsonschema submodule"
		git submodule update --init jsonschema
	else
		echo "Update jsonschema submodule"
		cd jsonschema
		git pull origin master
		cd ..
	fi
	#Libxsmm
	if [ ! -d libxsmm ]; then
		mkdir libxsmm
	fi
	if [ ! -f libxsmm/.git ]; then
		echo "Initialize libxsmm submodule"
		git submodule update --init libxsmm
		#Clean documentation to save space
		cd libxsmm
		rm -rf samples/
		rm -rf documentation/
		cd ..
	else
		echo "Update libxsmm submodule"
		cd libxsmm
		LOCAL=$(git rev-parse master)
		REMOTE=$(git rev-parse origin/master)
		if [ $LOCAL = $REMOTE ]; then
			echo "Up-to-date"
		else
			REBUILD_LIBXSMM=true
			git stash -q     #silently stash the changes (deleted directories)
			git pull origin master
			git stash pop -q #silently unstash the changes (deleted directories)
			rm -rf samples/       #delete potential new stuff
			rm -rf documentation/ #delete potential new stuff
		fi
		cd ..
	fi
	
	# build libxsmm
	if [ ! -d libxsmm/bin ] || [ ! -e libxsmm/bin/libxsmm_gemm_generator ]; then
		REBUILD_LIBXSMM=true
	fi
	if [ "$REBUILD_LIBXSMM" = true ]; then
		echo "Build libxsmm gemm generator"
		cd libxsmm
		make realclean
		make generator
		cd ..
	fi
else
	while getopts hsw opt; do
	case $opt in
		h)  echo "-h prints this message"
			echo "-s set submodules url to ssh"
			echo "-v set submodules url to https"
			exit -1;;
		s)  git config submodule.Submodules/Peano.url       git@gitlab.lrz.de:gi26det/Peano.git
			git config submodule.Submodules/jinja.url       git@github.com:pallets/jinja.git
			git config submodule.Submodules/markupsafe.url  git@github.com:pallets/markupsafe.git
			git config submodule.Submodules/attrs.url       git@github.com:python-attrs/attrs.git
			git config submodule.Submodules/pyrsistent.url  git@github.com:tobgu/pyrsistent.git
			git config submodule.Submodules/jsonschema.url  git@github.com:Julian/jsonschema.git
			git config submodule.Submodules/libxsmm.url     git@github.com:hfp/libxsmm.git
		w)  git config submodule.Submodules/Peano.url       https://gitlab.lrz.de/gi26det/Peano.git
			git config submodule.Submodules/jinja.url       https://github.com/pallets/jinja.git
			git config submodule.Submodules/markupsafe.url  https://github.com/pallets/markupsafe.git
			git config submodule.Submodules/attrs.url       https://github.com/python-attrs/attrs.git
			git config submodule.Submodules/pyrsistent.url  https://github.com/tobgu/pyrsistent.git
			git config submodule.Submodules/jsonschema.url  https://github.com/Julian/jsonschema.git
			git config submodule.Submodules/libxsmm.url     https://github.com/hfp/libxsmm.git
	esac
	done
fi

# move back to where the script was called
cd "$currentLocation"




