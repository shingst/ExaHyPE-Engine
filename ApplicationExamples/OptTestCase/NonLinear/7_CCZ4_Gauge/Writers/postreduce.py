#!/usr/bin/env python
# written 2017-03-14 by SvenK for ExaHyPE.

"""
Reduce TimeSeriesReductions written by several MPI ranks in a post processing step.

Usage: python2 postreduce.py ../path/to/files-rank*.asc
  or:  cd ../path/to; ./postreduce.py

This is a Python2 + numpy script which makes an educated guess which files to
merge together. It uses multiprocessing to exploit all cores for speed up
slurping the ASCII files (map phase). The reduction is then done by definitions
here in the Python script -- make sure they are similar to what you have defined
in your TimeSeriesReductions.h for the C++ part.
"""

import numpy as np

# all these are python standard builtins
from multiprocessing import Pool
from glob import glob
from collections import defaultdict, namedtuple
from sys import argv, exit
from os import path
import re
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

# This comes directly from TimeSeriesReductions.h.
# Note that the column ordering is read off the first input file
# and is thus arbitary in columns.
columns = {}
Col = namedtuple('Col', ['format', 'reduction'])

# handy reductions
noreduction = None
sums = lambda c: lambda r: np.sum(r[c], axis=0)

columns['plotindex'] = Col('%.0f', noreduction)
columns['time']      = Col('%e', noreduction)
columns['l1norm']    = Col('%e', sums('l1norm'))
columns['l2norm']    = Col('%e', sums('l2norm'))
columns['min']       = Col('%e', lambda r: np.min(r['min'], axis=0))
columns['max']       = Col('%e', lambda r: np.max(r['max'], axis=0))
columns['avg']       = Col('%e', lambda r: np.average(r['avg'], weights=r['numelements'], axis=0))
columns['numelements'] = Col('%.0f', sums('numelements'))
columns['numnan']    = Col('%.0f', sums('numnan'))

# columns which must be the same in all ranks
same_columns = { 'plotindex', 'time' }

# read argv for beginning
allfiles = argv[1:]

if not allfiles:
	# read the working directory
	allfiles = glob("*-rank-*.asc")
	logger.info("Read %i reductions ascii files from the current working directory" % len(allfiles))

if not allfiles:
	logger.error("No files found.")
	logger.info("Usage:")
	logger.info(__doc__)

# search for fields belonging together
fields = defaultdict(list) # mapping field name => list of input files
for fname in allfiles:
	if not path.isfile(fname):
		raise ValueError("File '%s' is not readable" % fname)
	parts = re.match(r"^(?P<fieldname>.+)-rank-(\d+)(?P<appendix>.asc)$", fname)
	assert parts, "Cannot read of rank number from filename '%s'" % fname
	
	# inserting this allows easier globbing on the command line
	globalfname = parts.group('fieldname') + "-global" + parts.group('appendix')
	fields[globalfname].append(fname)

def readasc(fname):
	# changes in the reductions CSV files go here
	return np.genfromtxt(fname, names=True)

def writeasc(fname, recarray, delimiter="\t", fmt='%.6e'):
	header = delimiter.join(recarray.dtype.names)

	# numpy.savetxt, at least as of numpy 1.6.2, writes bytes
	# to file, which doesn't work with a file open in text mode.  To
	# work around this deficiency, open the file in binary mode, and
	# write out the header as bytes.
	with open(fname, 'wb') as f:
		f.write(header + "\n")
		np.savetxt(f, recarray, fmt=fmt, delimiter=delimiter)

pool = Pool()

for fieldname, rankfilenames in fields.iteritems():
	ranktables = pool.map(readasc, rankfilenames)

	# The nested rec array is crucial for the reduction definitions
	ranktables = np.array(ranktables)

	assert len(ranktables) > 0, "Must have at least one rank to work with"
	globaltable = np.zeros_like(ranktables[0])
	foundcolumns = globaltable.dtype.names
	
	# assert the same times and number of elements in all ranks
	for colname in same_columns:
		assert np.all(ranktables[:][colname][0] == ranktables[:][colname]), "Inconsistencies in column '%s'" % colname
	for colname in foundcolumns:
		assert colname in columns, "I don't understand column '%s' as found in the input files" % colname

	# populate the initial values to each time row
	for colname in same_columns:
		globaltable[colname] = ranktables[0][colname]
	for colname, column in columns.iteritems():
		if column.reduction:
			globaltable[colname] = column.reduction(ranktables)

	formats = [columns[colname].format for colname in foundcolumns]
	writeasc(fieldname, globaltable, fmt=formats)
	
	logger.info("Reduced %i ranks to %s" % (len(rankfilenames), fieldname))


# Improvement ideas:
#
# 1. Don't overstress assertions, use a function wrapping logger.error + exit(1).
# 2. Use Engine/Misc/Postprocessing/lscompact.py to summarize the MPI ranks which
#    entered the global ASCII file in a one line comment like
#    "contains error-phi-rank-[0,52].asc, joined at 2017-03-14T23:00:00 at $(hostname)"
#    This is very helpful to detect noncontinous lists of files which otherwise could be 
#    a frequent mistake.
#
