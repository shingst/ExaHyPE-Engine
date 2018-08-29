#!/bin/bash
#
# Outsourced part of the exa.sh: Updating the repository.
#

SCRIPT="$(readlink -f $0)" # absolute path to update-installation.sh
ME="$(basename "$SCRIPT")" # just my name (update-installation.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
EXA="$SCRIPTDIR/exa.sh"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

verbose() { info $@; $@; } # only used for debugging
info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
subreq() { $SCRIPT $@; } # subrequest: Query another command for output
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory

cdroot;

set -e # <- this is important during the script execution!

# no arguments provided
[ $# -eq 0 ] && subreq help

# TODO: The idea was that one can call this script with several parameters
# and each is then done in order. However, what happens instead is that only
# the first argument is invoked.
# Thus the help message is misleading

for cmd in "$@"; do
	case $cmd in
		"exahype"|"bootstrap"|"all") # Does git pull, update-peano and update-toolkit
			info "Performing parallel quick standard exahype update"
			git pull
			subreq peano &
			subreq clean-toolkit & 
			# collecting the success of the two commands is not guaranteed.
			wait
			# make sure the Peano update succeeded. Sometimes when connections fail, the repo
			# remains in a bad situation and uses have to call "svn cleanup".
#			{ cd Peano && svn info; } || fail "Peano is assumably locked and the update FAILED. You have to manually clean it up."
			# We don't do libxsmm as it's not needed for most applications
			info "Finished updating ExahyPE (+peano, +toolkit) in $GITROOT"
			;;
		"peano") # Updates the Peano subversion repository
			info "Updating Peano"
			cd Peano
			exec ./updatePeano.sh
			;;
		"toolkit") # Compiles the toolkit with ant and javac
			info "Creating Toolkit"
			cd Toolkit
			exec ./build.sh
			;;
		"clean-toolkit") # Workaround for the toolkit build
			info "Creating a clean toolkit"
			cd Toolkit
			rm -rf src/eu/exahype/node*
			exec ./build.sh
			;;
		"toolkit-jar") # Just rebuild jar with new template files in toolkit
			info "Creating Toolkit JAR"
			cd Toolkit/src
			make dist
			;;
		"toolkit-fastmake") # no clean, just recompile and make jar
			info "Toolkit make and dist"
			cd Toolkit/src
			make compile && make dist
			;;
		"libxsmm") # Checkout or recompile libxsmm
			info "Cloning or updating Libxsmm"
			[[ -e Libxsmm ]] || git clone https://github.com/hfp/libxsmm.git Libxsmm
			# libxsmm is a large repository, instead of cloning, downloading
			# https://github.com/hfp/libxsmm/archive/master.zip
			# would be an option
			cd Libxsmm
			git pull
			make clean
			exec make generator -j4
			;;
		"full-eulerflow") # Update exahype, clean, recompile and run the EulerFlow application
			info "Making a full update and cleaning of the standard EulerFlow application"
			# this will be no more neccessary once we have out-of-tree builds.
			StandardApp="EulerFlow"
			subreq exahype
			cd $($EXA find-appdir $StandardApp)
			make clean || true # don't break if never run before
			$EXA compile-run $StandardApp
			;;
		""|"help"|"--help")
			info "Allows updating parts of this repository"
			info "Valid commands are: exahype/bootstrap/all, peano, toolkit, clean-toolkit, toolkit-jar, toolkit-fastmake, libxsmm"
			info "and all of it's combinations, whitespace seperated"
			;;
		*)
			fail "Could not understand command '$cmd'"
			;;
	esac
done
