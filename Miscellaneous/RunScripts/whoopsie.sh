#!/bin/bash
#
## This is a more generic version of Misc/BuildScripts/whoopsie-paster.sh, which is supposed
## to dump all kind of debug-related data. Instead, this script tries to do so in terms 
## of ExaHyPE executables.
##
## The short idea of this script is to dump all kind of information in a text file and upload
## it to a pastebin service such as sprunge.us or termbin.com. If this is not possible, instead
## provide the user the generated file for dumping on his own.
##
## Usage examples:
##
##    ../path/to/whoopsie.sh  ./ExaHyPE-Euler ../EulerFlow.exahype
##    exa whoopsie ./ExaHyPE-Euler ../EulerFlow.exahype
##    exa whoopsie mpirun -np 10 ./ExaHyPE-Euler ../EulerFlow.exahype
##
## Minimalistic demonstrators for showing how success and failure is captured:
##
##    exa whoopsie true
##    exa whoopsie false
##
## Written by SvenK at 2017-11-06.
#

runlog="run.log"    # where (only) the executable output goes
dumplog="dump.log"  # where also the whoopsie analysis stuff goes
# wipe the files
for f in $runlog $dumplog; do echo -n > $f; done;

# self documenting: if no parameter given
if [ $# -eq 0 ]; then grep -E '^##' "$0"; exit; fi

have() { which $@ 2>/dev/null >/dev/null; } # checks on the path, not bash builtins

runscripts="$(dirname "$0")"
buildscripts="${runscripts}/../BuildScripts"
exa="${buildscripts}/exa.sh"

# Turn off buffering in pipes. Set this to the empty string if your
# system does not have the command "stdbuf"
if have stdbuf; then unbuf="stdbuf -i0 -o0 -e0"; else unbuf=""; fi

# Use like "something | $teelog": Dumps both into all logfiles and stdout.
teelog="$unbuf tee -a $dumplog $runlog"

# output
hiddenexec() { $@ >> $dumplog 2>&1; } # execute something and put all output into dumplog
log() { hiddenexec echo $@; } # log something

# semantics of our funny minimalistic markup language
markup() { log "#{$1} ${@:2}"; }
heading() { markup "HEAD" $@; }
logcmd() { markup "CMD" $@; }
highlight() { markup "HIGHLIGHT" $@; }
logitem() { markup "ITEM" $@; } # list item
spacing() { log; } # just some lines of spacing
beginbracket() { log "#{BEGIN} $@"; }
endbracket() { log "#{END} $@"; }
bracketexec() { beginbracket $@; hiddenexec $@; endbracket $@; }

# composita
verbose() { logcmd $@; bracketexec $@; } # verbose a command and execute it
log_file() { bracketexec cat $@; } # log a text file

heading "Welcome to the ExaHyPE whoopsie runtime problem catcher output"
log "I accompany the user $(whoami) on the computer node $(hostname) at"
log "time $(date) by running ExaHyPE. Running ExaHyPE is not easy."
spacing
log "So the user starts with an environment given by"
spacing
verbose $exa check
spacing
heading "Running the actual command"

logcmd $@
# Try to run the commands, capture also the output. If it finishes: fine.
ulimit -c unlimited  # ulimit is bash builtin
set -o pipefail  # needed for detecting failure in pipes
beginbracket $@
if 2>&1 $@ | $teelog; then
	log "Finished successfully."
	exit
fi
endbracket $@

# else: Log all possible stuff
echo "Whoopsie: Catching a failed command to $PWD/$dumplog"
highlight "The command FAILED. " # with return value $? <= not the return value of the program. Needs Pipe inspection.
spacing
log "In the following, various environmental and command-related output will be collected."
spacing
heading "Post mortem command inspection"
spacing

# collect information about the run:
is_mpi="false"
specfile_path="unknown"

for part in $@; do
	if [[ $part == *"mpi"* ]] && have which && which $part &>/dev/null; then
		logitem "$part is something like mpirun or mpiexec. When I also find the specfile, I can run the Peano domaincomposition on the output."
		spacing
		log "This is where $part resolves to:"
		verbose which $part
		is_mpi="true"
	elif [[ -e $part ]]; then
		if [[ -x "$part" ]] && ./$part --help 2>&1 | grep -qi exahype; then
			logitem "$part is an ExaHyPE executable. This is how it was compiled:"
			spacing
			verbose $part --version
		elif have file; then
			if file $part | grep -qi ascii; then
				logitem "$part is a text file. Here are it's contents:"
				log_file $part
				if grep "exahype-project" $part && grep "computational-domain"; then
					logitem "$part is also very likely an exahype specfile. I will remember it for later use."
					specfile_path="$part"
				fi
			else
				log "$part is a file but I don't know what's inside. This is what I can learn about it:"
				verbose file $part
			fi
		else
			logitem "$part is a file but I miss the Unix tool 'file' to look into it."
		fi
	elif have which; then
		# only check parts which don't go like "-foo" or "--foo" because which interpretes this as parameters
		if ! [[ $part == -* ]] && which $part &>/dev/null; then
			logitem "$part is on the PATH. It resolves to:"
			verbose which $part
		else
			logitem "'$part' does not seem to be on the PATH nor it is a file. Probably a plain old parameter."
		fi
	else
		logitem "'$part' is not a file and I have no tools (no 'which') to find out what it is."
	fi
done

spacing
heading "MPI post mortem performance analysis"
if [[ $is_mpi == "true" ]]; then
	log "I detected an MPI call in this run.";
	if ! [[ $specfile_path == "unknown" ]] && [[ -e $specfile_path ]]; then
		log "And I found the specfile at $specfile_path"
		dimension=$(grep -m1 -A3 'computational-domain' $specfile_path | grep 'dimension' | cut -d'=' -f2 )
		domainsize=$(grep -m1 -A3 'computational-domain' $specfile_path | grep 'width' | cut -d'=' -f2 | tr -d ',')
		domainoffset=$(grep -m1 -A3 'computational-domain' $specfile_path | grep 'offset' | cut -d'=' -f2 | tr -d ',')
		$exa peano-analysis -dimension $dimension -domainoffset $domainoffset -domainsize $domainsize $runlog
		# TODO: collect output to files.
		highlight "I run the peano-analysis -> collecting output -> todo"
	else
		highlight "However, I could not find the specfile (it is: $specfile_path)."
	fi
else
	log "Did not detect an MPI call in this run."
fi

spacing
spacing
heading "The current environment: "
spacing
verbose env
spacing
heading "Inspection of local changes in the installation:"
verbose $buildscripts/installation-status.sh

log
if [[ -e make.log ]]; then
	heading "I have found a make.log which probably explains how this was built:"
	log_file make.log
else
	log "No make log found. Are there other log files? Let's check:"
	verbose ls *.log
fi

log "$0 is finished."

uploadWith() {
	echo "Uploading the file ${dumplog}, please send the following link to your collaborators:"
	cat $dumplog | $@ || { echo "Failure with uploading! Please just send the file $PWD/${extendedlog} as an email attachment to your collaborators."; }
}

if have curl; then
	uploadWith 'curl -s -F sprunge=<- http://sprunge.us'
elif have nc; then
	uploadWith 'nc termbin.com 9999'
else
	echo "Please send the file $PWD/${dumplog} as an email attachment to your collaborators."
fi

