#!/bin/bash
#
# This is some helper script to paste the make.log and additional information
# to some paste service.
#
# Several options:
#    <command> | curl -F 'sprunge=<-' http://sprunge.us
#    <command> 2>&1 | curl -F 'f:1=<-' ix.io
#    echo just testing!  | nc termbin.com 9999
#    dmesg | curl -F "upfile=@-" https://clitxt.com
# all returning the generated URL.
#
# A more structured pastebin option are https://gist.github.com/ which can
# easily be used anonymously but with some JSON work. Advantage: Post several
# files in one git commit.
# See API at https://developer.github.com/v3/gists/#create-a-gist
# Usage from command line with the given JSON payload stored in test.json is
#    curl -X POST -d @test.json https://api.github.com/gists  --header "Content-Type:application/json"
# but requires json output parsing (eg. with jq). It's much easier to use
# eg. python here to do the submission.
#
# For the time begin, this shell script uses the plain text pastebins only.
# Feel free to extend for github if you feel the need to do so.
#

# TODO: Improve documentation, consolidate with ../RunScripts/whoopsie.sh and ./analyze-git.sh

extendedlog="extended.log"

# wipe the file
echo -n > $extendedlog

log() { echo -e $@ >> $extendedlog; }
have() { which $@ 2>/dev/null >/dev/null; } 

log "This is an ExaHyPE 'whoopsie' compiler post mortem log."
log "It is suitable for showing the errors to other group members."
log "Read more about ExaHyPE at http://www.exahype.eu"
log
log "Compilation failed for $(whoami) at the computer named $(hostname)"
log "at $(date)."
log "We are in the directory $(pwd)"
log "and these are our main ExaHyPE variables:"
log
log "COMPILER=$COMPILER"
log "SHAREDMEM=$SHAREDMEM"
log "DISTRIBUTEDMEM=$DISTRIBUTEDMEM"
log "MODE=$MODE"
log "TBB_INC=$TBB_INC"
log "TBB_SHLIB=$TBB_SHLIB"
log
log "These are variables related to Svens exa compiler suite:"
log
log "SKIP_TOOLKIT=$SKIP_TOOLKIT"
log "CLEAN=$CLEAN"
log
log "We are basically at this commit:"
log "$(git show --oneline -s)"
log
log "And all modifications in current directory are:"
git status -s . >> $extendedlog

# This may expose private information and also adds lots of rubbish:
#log
#log "Actually, if you really want, here is all the environment I am using:"
#log "$(env)"
log
log "So, now let's get into something more interesting. Here is the make.log:"
log
cat make.log >> $extendedlog
log
log "Bye."

uploadWith() {
	echo "Uploading the file ${extendedlog}, please send the following link to your collaborators:"
	cat $extendedlog | $@ || { echo "Failure with uploading! Please just send the file $PWD/${extendedlog} as an email attachment to your collaborators."; }
}

if have curl; then
	# sprunge seems to have stopped working
	# uploadWith 'curl -s -F sprunge=<- http://sprunge.us'
	uploadWith 'curl -s -F clbin=<- https://clbin.com'
elif have nc; then
	uploadWith 'nc termbin.com 9999'
else
	echo "Please send the file $PWD/${extendedlog} as an email attachment to your collaborators."
fi

