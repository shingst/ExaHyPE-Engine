#!/bin/bash

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

git clone git@bitbucket.org:svek/mexa.git .

echo "Feeding in relevant commits from $ExaHyPE_Engine"

git remote add -f exahype-engine $ExaHyPE_Engine

git filter-branch --prune-empty --subdirectory-filter  Miscellaneous/MetaSpecfile/  exahype-engine/master

msg="Merge changes made in ExaHyPE-Engine repository into the standalone repository of mexa."

echo "Merging..."

git merge -X theirs -m"$msg" exahype-engine/master

echo "Pushing..."

git push

echo "Done. See $PWD for what happened (maybe: git log)."
