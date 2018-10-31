#!/bin/bash
#
# A wrapper around typing "make"  inside an ExaHyPE application
# directory, featuring
#
# 1. environment parameter driven behaviour (just like ExaHyPE's Makefile)
# 2. Invocation of the toolkit
# 3. Patching of faulty files where the toolkit is broken
# 4. Cleaning (also partially) before make
# 5. Logging the output of make, also passing to whoopsie in case of failure
# 6. Invoking auto parallelized make.
#
# To use, go to the ExaHyPE app you want to compile, and instead of typing
#   > make -j4 2&>1 | tee make.log
# Just type
#   > ../path/to/buildScripts/compile.sh
# Or just use exa:
#   > exa compile YourApp
# from anywhere.
#
# (c) 2016 ExaHyPE, Sven K

buildscripts="$(dirname "$0")"

fail() { >&2 echo $@: $@; exit -1; }
verbose() { >&2 echo $@; $@; }
has() { type $@ &>/dev/null; } # a way to check if command is available

# if SPECFILE is unset   AND  no cmd line args passed
if [ -z ${SPECFILE+x} ]  &&   [ $# -eq 0 ]; then
	# try to guess the specfile location by reverse guess
	SPECFILE=$($buildscripts/exa reverse-spec "$(PWD)")
fi

# options for the Make systems
DEFAULT_COMPILER="GNU"
DEFAULT_SHAREDMEM="None"
DEFAULT_DISTRIBUTEDMEM="None"
DEFAULT_MODE="Asserts"

# options controlling how this script works
DEFAULT_CLEAN="None"
DEFAULT_SKIP_TOOLKIT="No"
DEFAULT_MAKE_NPROC="$(has nproc && nproc || echo 1)"

# Allow to store information how to compile on specific machines
source $buildscripts/load-clusterconfig.sh

# all default variables can be overwritten by specifying them as
# environment variables

# if one of these two variables is unset:
if [ -z ${SPECFILE+x} ] || [ -z ${ABSCODEDIR+x} ]; then
	>&2 echo "Please specify SPECFILE and ABSCODEDIR or at LEAST ABSCODEDIR so I know where I shall run."
	exit -1
fi

CLEAN=${CLEAN:=$DEFAULT_CLEAN}
SKIP_TOOLKIT=${SKIP_TOOLKIT:=$DEFAULT_SKIP_TOOLKIT}
MAKE_NPROC=${MAKE_NPROC:=$DEFAULT_MAKE_NPROC}

# info: if you run into trouble here, move the echo statements from below up here.

# go to ExaHyPE-Engine root directory (used to be Code/ in former days)
verbose cd "$ABSCODEDIR" || fail "Cannot compile as there is no ABSCODEDIR=$ABSCODEDIR"
[[ -e "$SPECFILE" ]] || fail "Cannot find specfile $SPECFILE in $PWD";
PROJECTNAME=$(grep '^exahype-project' "$SPECFILE" | awk '{ print $2; }')
APPDIR=$(cat "$SPECFILE" | awk 'BEGIN{r=1} /output-directory/{ r=0; print $4; } END{ exit r}' || fail "Failed to determine Application output directory in SPECFILE=$SPECFILE" )

# Logging all further invocations of the toolkit, etc. to make.log
unbuf="stdbuf -i0 -o0 -e0" # turn off buffering in pipe
exec &> >($unbuf tee "$APPDIR/make.log")

echo "$0 running with"
echo " SPECFILE = $SPECFILE"
echo " ABSCODEDIR = $ABSCODEDIR"
echo " CLEAN = $CLEAN"
echo " CLUSTERNAME = ${CLUSTERNAME:=-not set-}"
echo " APPDIR = $APPDIR"
echo " PROJECTNAME = $PROJECTNAME"
echo " SKIP_TOOLKIT = $SKIP_TOOLKIT"

export COMPILER=${COMPILER:=$DEFAULT_COMPILER}
export SHAREDMEM=${SHAREDMEM:=$DEFAULT_SHAREDMEM}
export DISTRIBUTEDMEM=${DISTRIBUTEDMEM:=$DEFAULT_DISTRIBUTEDMEM}
export MODE=${MODE:=$DEFAULT_MODE}

# you can amend on this
#  export TBB_INC=/usr/include/tbb
#  MPI_LDFLAGS="$(mpicc -showme:link)" # eigentlich: mpiicpc -showme:link ...
#  export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

#echo -e " COMPILER=$COMPILER"
#echo -e " SHAREDMEM=$SHAREDMEM"
#echo -e " DISTRIBUTEDMEM=$DISTRIBUTEDMEM"
#echo -e " MODE=$MODE"

echo -e "at $(date) on $(hostname) as $(whoami)"
echo -e

set -e

# run the toolkit on this application
if [[ $SKIP_TOOLKIT == "Yes" ]]; then
	echo -e "Skipping toolkit invocation as requested";
else
	echo -e "Running toolkit"

	# todo: 
	#echo -e "Working around defect Makefiles etc"
	#rm $APPDIR/Makefile
	#could also delete KernelCalls.cpp, $APPNAME_generated.cpp, etc.

#	verbose java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $SPECFILE || { >&2 echo "Failure when running the toolkit"; exit -1; }
	verbose Toolkit/toolkit.sh $SPECFILE || { >&2 echo "Failure when running the toolkit"; exit -1; }
fi

cd -

# Allow to store information how to compile a specific application
APP_INFO_FILE="projectpaths.cfg"
if [[ -e "$APP_INFO_FILE" ]]; then
	echo "Loading application specific configuration from $APP_INFO_FILE"
	source $APP_INFO_FILE
	# the app config file is supposed to set stuff like:
	export PROJECT_CFLAGS="${PROJECT_CFLAGS:=}"
	export PROJECT_LFLAGS="${PROJECT_LFLAGS:=}"
	export PROJECT_LINK="${PROJECT_LINK:=}"
	echo " PROJECT_CFLAGS = $PROJECT_CFLAGS"
	echo " PROJECT_LFLAGS = $PROJECT_LFLAGS"
	echo " PROJECT_LINK = $PROJECT_LINK"
	
else
	echo "No application specific configuration ($APP_INFO_FILE) found."
fi

# plausability check
[[ -e Makefile ]] || { echo -e "Could not find Makefile in $PWD. Probably the toolkit failed."; exit -1; }

case $CLEAN in
	"None") echo -e "No cleaning before building."
		;;
	"Clean") 
		verbose make clean
		;;
	"Lightweight") echo -e "Lightweight clean"
		# find also object files in subdirectories
		verbose find . -iname '.o' -exec rm {} \;
		# and also cleanup Fortran modules
		verbose find . -iname '.mod' -exec rm {} \;
		;;
esac

# Workaround the broken makefile system
echo -e "Fixing Makefile after toolkit run"
sed -i "s/SHAREDMEM=.*/SHAREDMEM=$SHAREDMEM/" Makefile
sed -i "s/DISTRIBUTEDMEM=.*/DISTRIBUTEDMEM=$DISTRIBUTEDMEM/" Makefile

# Workaround for files overwritten by toolkit
for patchfile in "${PROJECTNAME}_generated.cpp"; do # add files seperated by whitespace as in "foo" "bar"
	if [ -e $patchfile ]; then
		if git ls-files $patchfile --error-unmatch &>/dev/null; then
			echo -e "Patching $patchfile with committed version to overwrite toolkit changes"
			git checkout $patchfile
		else
			echo -e "Patching: File $patchfile not under version control"
		fi
	else
		echo -e "Don't patch $patchfile as not there"
	fi
done

# Workaround for broken dependency build system:
# Generate well known fortran modules if present

# This works, but cannot access the EXAHYPE_CFLAGS or similar. Instead, I made
# a small change in the Makesystem so Modules compile first. Probably.

#for fmodule in Parameters.f90 typesDef.f90; do
#	if [ -e $fmodule ]; then
#		echo -e "Precompiling $fmodule as otherwise build fails"
#		FORTFLAGS="-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none"
#		verbose gfortran $FORTFLAGS -c $fmodule
#	fi
#done

# Force GCC to colour output even in the redirected output.
if [[ $COMPILER == "GNU" ]]; then
	gcc_force_color_output="-fdiagnostics-color=always"
	export PROJECT_CFLAGS="${PROJECT_CFLAGS} $gcc_force_color_output"
fi

set -o pipefail # fail if make fails

verbose make -j $MAKE_NPROC || {
	echo -e "Make failed!";
	$buildscripts/whoopsie-paster.sh
	exit -1;
}

echo -e "Compile.sh finished successfully"

