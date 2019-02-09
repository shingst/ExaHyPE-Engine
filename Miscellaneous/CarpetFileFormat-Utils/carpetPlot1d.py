#!/usr/bin/env python
from pylab import np, plt
import h5py, sys, os, re

"""
This is a minimal demonstrator program to show how to read
and plot Carpet HDF5 1d data in Python.

Usage:
  ./carpetPlot1d.py  yourStuff.x.h5  anyOtherStuff-*.h5 ...
"""

within = lambda key, lst: [ k for k in lst if key in k ]

def x_coordinates_for(dataset):
	origin = np.asscalar(dataset.attrs["origin"])
	delta = np.asscalar(dataset.attrs["delta"])
	N = len(dataset)
	return origin + np.arange(N) * delta

for fname in within("h5", sys.argv):
	with h5py.File(fname) as hf:
		for dataset in hf["Parameters and Global Attributes/Datasets"]:
			keys = within(dataset, hf)
			timesteps = np.unique([hf[k].attrs["timestep"] for k in keys])
			for ts in timesteps:
				fout = "%s-field-%s-it%d.png" % (os.path.splitext(fname)[0], re.sub("::", "-", dataset), ts)
				print "Plotting %s..." % fout
				data = [hf[k] for k in keys if hf[k].attrs["timestep"] == ts]
				x = np.concatenate(map(x_coordinates_for, data))
				y = np.concatenate(data)

				plt.clf()
				plt.title("%s - timestep %d" % (dataset, ts))
				plt.plot(x,y, "-o")
				plt.savefig(fout)
