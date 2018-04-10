#!/bin/bash
#
# exa is a versatile CLI entry point for various scripts.
# Install it with a softlink from your /usr/local/bin or ~/bin to
# make use of the full delocalized power.
#
#
# Idee: Ebenfalls als bashrc script nicht nur fuer completion sondern auch slimmeres gefuehl
#       beim Tabben und und "exa cd" oder "exa paths" erlaubt setzen von umgebungsvariablen direkt
#       drinnen.
#       Alternativ: exa shell mit Parametern und veraenderter bashrc? finde ich weniger sinnvoll.
#
# (c) 2016 ExaHyPE - by SvenK


# Compile vs build
# ================
#
# A note about the "compile" vs "build" idiom: All compile-* commands are
# shiny wrappers around compile.sh and thus around make, using ExaHyPE's
# Makefile.
#
# In contrast, all "build-*" commands are wrappers around the out of tree
# system which allows to build basically in copies of the repository and
# transfering the executables back. Thus, we adopt there the language of
# "builds" as configurations of files which can be setup, synced, compiled
# or otherwise manipulated after an initial setup.
#

SCRIPT="$(readlink -f $0)" # absolute path to exa.sh
ME="$(basename "$SCRIPT")" # just my name (exa.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)

shopt -s expand_aliases
alias pop='set -- "${@:2}"'

pop # pop first parameter

err() { >&2 echo $@; }
verbose() { info $@; $@; } # only used for debugging
info () { err $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { err $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
subreq() { $SCRIPT $@; } # subrequest: Query another command for output
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory
# set $APPNAME or die
getappname() {
	# First, check for USER_CALLING_DIR and make some guess with the "find reverse" tool.

	APPNAME="$PAR";
	if [ -z "$APPNAME" ]; then
		# Appname not given as argument. Try instead to obtain it from PWD.
		### TODO: Allow a call to APPNAME="$(subreq find reverse-pwd)"
		###       but it's not yet working.
		abort "Usage: $0 $CMD <AppName>";
	fi
}
# To be called at the beginning of a command. Change to application directory and
# obtain the application name afterwards. *If* no application name was given, try
# to guess the application name from the current directory (pwd).
cdapp() {
	export USER_CALLING_DIR="$PWD"
	cdroot;
	getappname;
	cd $(subreq find appdir "$APPNAME") || abort "Could not go to app";
}

# some paths to the exa helper scripts
BuildScripts=$GITROOT/Miscellaneous/BuildScripts # == $SCRIPTDIR
Postprocesing=$GITROOT/Miscellaneous/Postprocessing

case $CMD in
	"update") # Update the repository or dependencies. Use "--help" for help.
		exec $BuildScripts/update-installation.sh $@
		;;
	"find") # Locate apps, spec files, compiled files, etc.
		exec $BuildScripts/find-app.sh $@
		;;
	"clusterconfig"|"config")  # Load cluster specific settings. Usage: "eval $(exa config)" or "exa config iboga-gcc-tbb"
		cdroot
		# of course there is not much purpose in sourcing this as exa.sh is currently
		# not be intended to be sourced. What we could do here is to echo the ENV
		# so it can be used like "source <(exa config)" or similar.
		source $BuildScripts/load-clusterconfig.sh $@
		# Currently, users can at least "eval $(exa config something)"
		# or just call "$(exa config something)" from their command line
		;;
	"toolkit") # Run the toolkit for an application, without compiling
		cdroot; getappname
		SPECFILE="$(subreq find specfile "$APPNAME")" || abort "Could not find specfile: $SPECFILE"
		info "Running ExaHyPE.jar on $SPECFILE"
		exec java -jar Toolkit/dist/ExaHyPE.jar --not-interactive $SPECFILE
		;;
	"compile") # Invokes the toolkit and compilation of an application
		cdapp;
		export SPECFILE="$(subreq find specfile "$APPNAME")" || abort "Could not retrieve specfile of $APPNAME, I just got $SPECFILE"
		export ABSCODEDIR="$GITROOT"
		$SCRIPTDIR/compile.sh "$SPECFILE"
		;;
	"compile-run") # A shorthand for compiling and running an application
		getappname
		subreq compile "$APPNAME"
		subreq run "$APPNAME"
		;;
	"compile-poly") # Compile for different polynomial orders (as basis for convergence studies).
		pop
		cdapp; $SCRIPTDIR/compile-for-polyorder.sh $@
		;;
	"compile-poly-all") # Compile for polynomial orders 2,3,4,5,6,7,8,9, serially
		cdroot; getappname
		# do NOT parallelize the loop as the build system does not allow
		# Use build-poly-all for a parallel version
		set -e
		for p in 2 3 4 5 6 6 7 8 9; do subreq compile-poly $APPNAME $p; done
		;;	
	"batch-compile") # build-compile with batch parameter support.
		exec $SCRIPTDIR/batch-compile.sh $@
		;;
	"make") # compile without invoking the toolkit
		cdapp
		export SKIP_TOOLKIT="Yes"
		export CLEAN="${CLEAN:=Lightweight}" # do no heavy cleaning
		subreq compile "$APPNAME"
		;;
	"build-init") # Initialize an out of tree build, ie. dont sync. Parameters: <AppName> <BuildName>
		cdroot; getappname; buildName=$2
		SPECFILE="$(subreq find specfile $APPNAME)" || abort "Could not find specfile: $SPECFILE"
		exec $SCRIPTDIR/setup-out-of-tree.sh $SPECFILE $buildName
		;;
	"build-setup") # Setup an out of tree build, ie. init and sync. Parameters: <AppName> <BuildName>
		cdroot; getappname
		eval $(subreq build-init $@) || abort "Could not init out of tree build."
		$oot_outdir/sync.sh
		;;
	"build-find") # Find the location of a build from parameter <BuildName>
		cdroot; buildName=$1
		[[ -e "Builds/build-$buildName" ]] || abort "Could not find build '$buildName'. List of builds: $(ls Builds/)"
		echo Builds/build-$buildName
		;;
	"build-sync") # Sync an out of tree build. Parameters: <BuildName>
		cdroot; buildName=$1
		buildloc=$(subreq build-find $buildName) || abort "Could determine build location"
		exec $buildloc/sync.sh
		;;
	"build-exec") # Execute a build executable. Parameters: <BuildName>
		# Can execute from anywhere. You might want to use as
		#   exa build-exec name-of-run path/to/some/specfile.exahype
		# to run in the PWD.
		buildName=$1
		buildloc=$(subreq root)/$(subreq build-find $buildName) || abort "Could determine build location"
		[[ -e $buildloc/oot.env ]] || abort "Cannot find build instance '$buildName'. Maybe execute 'exa build-compile $buildName' before?"
		source $buildloc/oot.env
		[[ -e $buildloc/$oot_binary ]] || abort "Cannot find build binary. Is the build finished?"
		exec $buildloc/$oot_binary ${@:2}
		;;
	"build-compile") # Setup and compile an out of tree build. Parameters: [AppName] <BuildName>
		cdroot; getappname; buildName=$2;
		set -e;
		subreq build-setup $APPNAME $buildName || abort "Could not setup build."
		buildloc=$(subreq build-find $buildName) || abort "Could determine build location"
		cd $buildloc && ./make.sh
		;;
	"build-poly") # Setup an oot build and compile for a given polynomial order. This is parallelizable.
		cdroot; getappname; pOrder=$2
		set -e
		[[ x$pOrder != "x" ]] || abort "Usage: <AppName> <pOrder>"
		buildName="p$pOrder"
		subreq build-setup $APPNAME $buildName || abort "Could not setup build for $buildName"
		buildloc=$(subreq build-find $buildName) || abort "Could determine build location"
		source $buildloc/oot.env
		cd $oot_builddir/$oot_appdir
		info "Compiling for p=$pOrder in $PWD"
		export CLEAN="Clean" # to avoid any side effects
		$SCRIPTDIR/compile-for-polyorder.sh $pOrder
		cdroot
		cp $oot_buildir/$oot_appdir/${oot_binary}-p$pOrder $oot_appdir/
		;;
	"make-clean") # clean an application
		cdapp
		make clean || fail "Cannot clean since toolkit did not run."
		;;
	"clean-at") # Clean all object files below some given path
		for f in $@; do
			if ! [[ -d $f ]] && [[ -e $f ]] && [[ -d $(dirname "$f") ]]; then
				echo "Cleaning at $(dirname $f)"
				subreq clean-at $(dirname "$f")
			else
				echo "Cleaning at $f"
				find $f -type f -iname \*.o -delete -print
			fi
		done;
		;;
	"cheat") # show the environment variables available for driving the build
		cdroot; cd $SCRIPTDIR;
		cat cheat-sheet.txt
		;;
	"check") # Tell which build flags we currently have in ENV
		err "ExaHyPE Makefile specific:"
		err "COMPILER:  ${COMPILER:=-not set-}"
		err "CC:        ${CC:=-not set-}"
		err "MODE:      ${MODE:=-not set-}"
		err "SHAREDMEM: ${SHAREDMEM:=-not set-}"
		err "DISTRIBUTEDMEM: ${DISTRIBUTEDMEM:=-not set -}"
		err "TBB_INC=${TBB_INC:=-not set -}"
		err "TBB_SHLIB=${TBB_SHLIB:=-not set -}"
		err
		err "Exa Project Makefile specific:"
		err "PROJECT_CFLAGS: ${PROJECT_CFLAGS:=-not set-}"
		err "PROJECT_LFLAGS: ${PROJECT_LFLAGS:=-not set-}"
		err
		err "Exa Build tool specific:"
		err "CLEAN:     ${CLEAN:=-not-set-}"
		err "SKIP_TOOLKIT: ${SKIP_TOOLKIT:=-not set-}"
		env | grep -iE '(exa|sim)' | grep -vE '^PWD|^OLDPWD'
		err
		err "Runtime specific:"
		err "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
		err
		err "Machine readable:"
		
		quoteenv() { perl -e 'foreach $k (sort(keys(%ENV))) { print "export $k=\"$ENV{$k}\"\n"; }'; }
		quoteenv | grep -iE '^export (COMPILER|CC|MODE|SHAREDMEM|DISTRIBUTEDMEM|CLEAN|SKIP_TOOLKIT|EXAHYPE_|PROJECT_)'
		;;
	"git") # passes commands to git
		cdroot; info "ExaHyPE Git Repository at $GITROOT"
		exec git $@
		;;
	"pwd") # give current working directory relative to ExaHyPE root	
		echo "Your directory:   $PWD"
		echo "ExaHyPE root dir: $GITROOT"
		fail "Calculation not yet implemented"
		;;
	"player") # Calls the 'exaplayer' plotting toolkit. Use "--help" for help.
		exec $Postprocessing/exaplayer.py $@
		;;
	"reader") # Calls the 'exareader' data conversion toolkit. Use "--help" for help.
		exec $Postprocessing/exareader.py $@
		;;
	"slicer") # Calls the 'exaslicer' VTK slicing toolkit. Use "--help" for help.
		exec $Postprocessing/exaslicer.py $@
		;;
	"peano-analysis") # Quickly start Peanos Domaincomposition analysis script.
		set -e
		exec python $GITROOT/Peano/peano/performanceanalysis/performanceanalysis.py $@
		;;
	"run") # quickly start an application inside it's directory. Cleans VTK files before.
		cdapp
		info "Starting $APPNAME with SKIP_TESTS." | tee run.log
		rm -f *.vtk *.log-file
		export EXAHYPE_SKIP_TESTS=TRUE
		ROOT=$(subreq root)
		BINARY=$(subreq find binary $APPNAME) && [[ -e "$ROOT/$BINARY" ]] || fail "Could not find binary ($BINARY)"
		SPECFILE=$(subreq find specfile $APPNAME)
		# the directory handling of this tool is really awkward.
		# Would prefer relative directories here.
		reducedbuf="stdbuf -oL -eL" # for quicker output, no 4k buffering
		$reducedbuf $ROOT/$BINARY $ROOT/$SPECFILE 2>&1 | $reducedbuf tee -a run.log
		;;
	"make-output-folders") # Create output folders for a given specfile
		cdapp
		ABS_SPECFILE=$GITROOT/$(subreq find specfile $APPNAME)
		echo "Creating folders in $(pwd)"
		exec $BuildScripts/../RunScripts/createOutputFoldersForSpecfile.sh $ABS_SPECFILE
		;;
	"sim") # lightweight simulation managament
		exec $BuildScripts/../RunScripts/sim.sh $@
		;;
	""|"help") # prints out help about the exa toolkit
		me=$(basename "$0")
		echo -e "$me: <command> [parameters]"
		echo -e "An ExaHyPE quick command helper."
		echo -e "It is operating with the ExaHyPE installation at $GITROOT"
		echo -e
		echo -e "Available commands:"
		echo -e
		cat $0 | grep -oE '^\s+"(.+)"\)' | tr -d '")|'
		echo -e
		echo -e "With their individual meanings:"
		echo -e
		# @TODO: Improve display of available formats
		cat $0 | grep -E '\)\s+#' | tr ')' ':' | tr -d '#' | column -c 2 
		echo -e
		;;
	"help-shortlist") # A list of commands, for shell completion
		cat $0 | grep -oE '^\s+"(.+)"\)' | tr -d '")|'"\n" | sed -e 's/[[:space:]]*/ /'
		;;
	"is") # prints out fortunes
		echo "cool"
		;;
	"todo") # prints out all files where TODO notes are inside
		cdroot; info "Probably all current files containing todo messages"
		find . -type f | grep -Ei '\.(cpp|C|h|f90|cc|tex)$' | xargs grep -i todo
		;;
	"root") # prints out the root of the ExaHyPE installation
		echo $GITROOT
		;;
	"whoopsie") # runs command, if fails, collects info and uploads to a pastebin
		exec $BuildScripts/../RunScripts/whoopsie.sh $@
		;;
	"mexa") # Runs the meta specfile file format utility
		exec $BuildScripts/../MetaSpecfile/mexa.py $@
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac

