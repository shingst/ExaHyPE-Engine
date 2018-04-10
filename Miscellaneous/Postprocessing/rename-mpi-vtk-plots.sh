#!/bin/bash
#
# This renames your files foo-7-rank-8.vtk to foo-rank-7-8.vtk.
# Usage:  scriptname foo
# where foo is the prefix of your files
#

rename -v "s,$1-([0-9]+)-rank-([0-9]+).vtk,$1-rank-"'$2-$1.vtk,' ${1}*.vtk
