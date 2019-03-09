#!/usr/bin/env python
#
# This script is part of the "Exa" Python Postprocessing tools.
# They are written in Python 2 with Python 3 compatibility in mind.
# Dependencies: VTK and numpy/scipy.
#
# This file is part of the ExaHyPE project.
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt

"""
exaslicer.py is both a standalone program as well a mini library to read VTK
ExaHyPE files, slice 3D down to 2D or 1D and output it as CSV or similar.

Important notes
===============

* This code does *not* interpolate. It really boils down to mask your
  selection criterion to the data available. If you select x==0 but there
  is data only on x=+-0.01, no interpolation will take place and you won't
  get any data.

* Therefore, this kind of "slicing" is only useful in uniform grids where you
  know exactly about where your data live. However, you can use the
  function close_mask to avoid at least looking up the coordinates.

Examples
========

# Just cut out the y-z plane where x has the value 0.5
# This gives you *no* data if there is no x value with 0.5.
x == 0.5

# However, inequalities will typically allow you more flexibility
y < 1e-5

# Use any numpy function, ie. for implementing tolerance
np.abs(x - 0.5) < 1e-10

# You can also combine masks (boolean operations)
# This gives you a 1D line in z direction:
(x == 0.5) & (y == 0.4)

# And this gives you ideally only one point
(x == 0.1) & (y == 1.0) & (z == 0.3)

# while this gives you probably two distinct planes
(x == 0.1) | (x == 0.9)

# You can also slice in time
(x == 0.5) & (t <= 10)

# Or you write for instance to get the z axis on any close coordinates,
# where close_mask is a function provided by the exa tools:
close_mask(x,0) & close_mask(y=0)

Notes
=====

* Slicing conditions are evaluated as python statements and can call any function.
  x,y,z,t are provided as numpy arrays, all with same shape (one dimensional lists)
  and it is expected that the result is a mask for similarl shape (i.e. same length).

* Slicing in time is not useful if you combine this with the ExaPlayer.

* Disable sorting the output in coordinates can help you to understand wether the
  slicing worked correctly.

* Use the 'coords' output format (output_coordinates) of the ExaWriter
  to get a quick view of the input and output coordinates.
  A typical workflow may be:

  $ exareader conserved-{1,2,3}.vtk --gridtype cells -o coords
  -> Read out some interesting coordinate, ie. x[114] = 0.525926
  $ exaslicer conserved-{1,2,3}.vtk --gridtype cells --sort-output --slice 'x==0.525926'

  in order to obtain a sliced and sorted CSV table.

"""

from __future__ import print_function
import sys, argparse, logging

import numpy as np
import matplotlib.pyplot as plt

from exahelpers import ExaFrontend, cleandoc, find_nearest, close_mask
from exareader import ExaReader, ExaWriter

logger = logging.getLogger("exaslicer")

# inspired by check1dsymmetry.py
timecoord = "time"
spacecoords = ["x", "y", "z"]
allcoords = ["time", "x", "y", "z"]

def slice(data, expr):
	# our approach is really straight forward: `expr` has to be correct python.
	t,x,y,z = data[timecoord], data['x'], data['y'], data['z']
	try:
		mask = eval(expr)
		if np.all(mask):
			logger.error("Nothing filtered, all elements survided")
			return None # trigger failure
		if not np.any(mask):
			logger.error("No element survived the filter.")
			return None
		return data[mask]
	except SyntaxError as e:
		logger.critical("Your expression is not valid Python.")
		logger.exception(e)
		return None		
	except NameError as e:
		logger.critical("Please only use t,x,y,z as variables to filter.")
		logger.exception(e)
		return None
	except Exception as e:
		logger.critical("There is something wrong with your expression.")
		logger.exception(e)
		return None


class Slicer:
	"""
	The ExaSlicer allows you to easily slice 3+1 dimensional data into pieces.
	This goes by filtering each data row with coordinate (t,x,y,z) according
	to your `slice` expression. These expressions are numpy masks.
	"""
	def __init__(self):
		pass
	
	def add_group(self, argparser):
		group = argparser.add_argument_group('slicing', description=cleandoc(self))

		group.add_argument('--slice', default=None, type=str,
		                   help="Slicing expression, e.g. 'x==0.001' or '(x==1e-3)&(y==5)'")
		group.add_argument('--no-sorting', action='store_true', default=True,
		                   help="Do not sort output in coordinates (t,x,y,z)")

	def apply_args(self, args, parser):
		# consistency check
		if not '==' in args.slice:
			logger.warn("Your expression '%s' is likely to be wrong. Right expressions look like 'x==0.5'.")

	@staticmethod
	def main(args=None):
		"""
		Just call this as Slicer().main(sys.argv) or omit the args
		"""
		
		frontend = ExaFrontend(epilog=__doc__,
			program_description='ExaHyPE data slicer')
		slicer = Slicer()#allow_slice_time=True)
		reader = ExaReader(files_as_main=True)
		writer = ExaWriter()
		frontend.add_module(reader)
		frontend.add_module(slicer)
		frontend.add_module(writer)

		args = frontend.parse_args(args)

		#args = reader.remove_prog_in_inputfiles_list(frontend.parser, args)
		sol = reader.read_files_as_requested(args)

		prev_shape = sol.shape
		logger.info("Slicing with mask solution = solution[%s]" % args.slice)	
		data = slice(sol, args.slice)
		if type(data) == type(None): # to avoid numpy casting to bool with error.
			frontend.parser.error("Slicing failed with your expression '%s'" % args.slice)

		logger.info("Sliced %s -> %s shaped array" % (prev_shape, data.shape))

		if not args.no_sorting:
			logger.info("Sorting output in order %s" % str(allcoords))
			data.sort(order=allcoords)

		writer.write_output_as_requested(data, args)

if __name__ == "__main__":
	Slicer.main() # sys.argv

