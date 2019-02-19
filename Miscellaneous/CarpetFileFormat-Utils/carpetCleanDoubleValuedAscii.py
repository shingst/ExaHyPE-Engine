#!/usr/bin/env python

import numpy as np, sys, re

"""
This is a small helper script to erase double valued data in ExaHyPE-generated
CarpetASCII output files. These data can look weird in case of when generated
by a FV solver, while DG solvers typically have better looking output.

Usage:
  ./carpetCleanDoubleValuedAscii.py  yourFile.x.asc AnotherFile.x.asc Another.xyz.asc
"""

def readFileIntro(fname):
	with open(fname) as fh:
		return "".join([ line for line in fh if re.match(r"^[^0-9]", line) ])
within = lambda key, lst: [ k for k in lst if key in k ]
def flatten(l): return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]
withoutComments = lambda fname: (line for line in open(fname,'r') if line[0] != '#')
nonempty = lambda lst: [ v for v in lst if v ]
toOrdinaryArray = lambda reca: reca.view('<f8').reshape(len(reca), len(reca.dtype))
recArrayMean = lambda reca: np.mean(toOrdinaryArray(reca), axis=0).view(reca.dtype)

for fname in within("asc", sys.argv):
	header = readFileIntro(fname)
	# if your output files have a name row, set names=True, else False
	table = np.genfromtxt(withoutComments(fname), names=True)
	for d in [dim for dim in "xyz" if dim in table.dtype.names ]:
		u, i = np.unique(table[d], return_inverse = True)
		double_rows = np.argwhere(np.bincount(i) > 1)
		delete_row_indices = []
		for rowi in double_rows:
			closeness_mask = np.isclose(table[d], table[rowi][d])
			candidates = np.argwhere(closeness_mask)
			if len(candidates) == 1: continue # apparently this is already handled
			survivor_index = candidates[0]
			delete_row_indices.append(list(candidates[1:].flatten())) # postpone to end
			table[survivor_index] = recArrayMean(table[closeness_mask])
		fixed_table = np.delete(table, np.unique(flatten(nonempty(delete_row_indices))))
		print "Merging %d double valued %s coordinates in %s..." % (len(delete_row_indices), d, fname)
		#assert len(np.unique(fixed_table[d])) == len(fixed_table[d]), "Uniqueness did not work"
		table = fixed_table # ready for next dimensional round

	outfname = fname + "-nonDoubleValued"
	with open(outfname, "w") as fh:
		print "Writing to %s" % outfname
		fh.write("# Double values sanatized by carpetCleanDoubleValuedAscii\n")
		fh.write(header)
		np.savetxt(fh, table)
	#np.savetxt(fname, table, header="# Double values sanatized by carpetCleanDoubleValuedAscii\n" + header)
	
