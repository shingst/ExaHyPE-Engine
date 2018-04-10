#!/usr/bin/env python
#
# This script enriches a given specfile with a combination of plotters.
# That can be used to compare the outcome of these plotters.
# It just prints the generated specfile which then can be further processed.
#
# Note that here we want each plotter to dump into its own output folder.
# ExaHyPE currently does not create folders when needed. Instead, use a
# bash onliner like
#
#   for output in $(grep vtk-output  ../GRMHD_cpp_Limiting3D_withAllPlotters.exahype | cut -d'=' -f2); do mkdir -p $(dirname "$output"); done;
#
# to create the neccessary output folders
#
# -- SvenK, 2017-11-29
#

from sys import argv
from itertools import product

mapping = "ConservedWriter" # the writer class used for comparing these plotters
base_specfile = argv[1] # use the filename as input argument
#base_specfile = "GRMHD_cpp_ADERDG-3D.exahype" # the base specfile to be overloaded

# all combinations of *k
combine = lambda *k: list(product(*k))
# flatten a potential 2d list, no flatten in tuples
def flatten(l):
    for el in l:
        if isinstance(el, list):# and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def plotterstring(plt, mapping=mapping):
	# plt is a tuple of a plotter configuration
	plotter_id = "::".join(plt)
	output_name = mapping + "-" + "-".join(plt)
	return """
		plot %(plotter_id)s %(mapping)s
			variables const = 23
			time      = 0.0
			repeat    = 0.0000001
			output    = ./vtk-output/%(output_name)s/%(output_name)s
		end plot
		""" % locals()

def print_inject_file(fname, magicline, injection):
	# python2 print, without newline.
	for line in open(base_specfile, "r"):
		if magicline in line:
			print injection,
		print line,

# all possible values:
formats = ("vtu","vtk","Peano","carpet")
basis = ("Cartesian", "Legendre")
storeposition = ("vertices", "cells")
layout = ("binary", "ascii", "hdf5")

plotters = list(flatten([
	combine(["vtk"],basis,storeposition,["binary","ascii"]),
	combine(["vtu"],basis,storeposition,["ascii"]),
	combine(["Peano"],basis,storeposition,["ascii","hdf5"]),
	combine(["vtk","vtu"],["Cartesian"],storeposition,["limited"],["binary","ascii"]),
	("carpet","cartesian","vertices","hdf5"),
	("vtk","patches","boxes","ascii"),
	("vtu","patches","boxes","ascii"),
]))

plotterstrings = map(plotterstring, plotters)

print_inject_file(base_specfile, "end solver", "\n".join(plotterstrings))
