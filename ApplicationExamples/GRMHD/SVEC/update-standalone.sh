#!/bin/bash
#
# This script tries to keep the standalone SVEC repository to be
# in sync with the code maintained at the Engine-ExaHyPE repository.
# To do so, it brutally merges all relevant commits into the
# standalone repository. These are filtered by the path name.
#
# Script based on Misc/MetaSpecfile/update-standalone.sh
# SvenK for ExaHyPE, 2018-02-15
#

set -e
cd "$(dirname $0)"

# ExaHyPE-Engine gitroot
ExaHyPE_Engine="$(git rev-parse --show-toplevel)"

if ! git remote -v | grep -q ExaHyPE-Engine.git; then
	echo "This script is supposed to Sync from ExaHyPE-Engine to standalone"
	exit 1
fi

tmpdir="$(mktemp -d)"
cd $tmpdir

echo "Checking out Standalone repository at $tmpdir"

git clone git@bitbucket.org:svek/svec.git .

echo "Feeding in relevant commits from $ExaHyPE_Engine"

git remote add -f exahype-engine $ExaHyPE_Engine

git filter-branch --prune-empty --subdirectory-filter  ApplicationExamples/GRMHD/SVEC  exahype-engine/master

msg="Merge changes made in ExaHyPE-Engine repository into the standalone repository of SVEC."

echo "Merging..."

git merge -X theirs -m"$msg" exahype-engine/master

echo "Pushing..."

git push

echo "Done. See $PWD for what happened (maybe: git log)."
