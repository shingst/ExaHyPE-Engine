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
A minimalistic interactive movie player for 1D or 2D data.
currently only displays x-y layer or x axis.

Example invocation on command line:

   ./exaplayer.py solution-*.vtk -i vtk -d 1 -s -I

This will interactively show up a window with a range slider
where you can look at a 1D slice over the data over the time.

You can do the same and look into the data in an interactive
IPython session. This allows you fine control eg. for creating
a movie, start with 'ipython -i' and do

>>> from exareader import exareader
>>> from exaplayer import Play1D
>>> from glob import glob
>>> import matplotlib.pyplot as plt
>>> fnames = glob('*.vtk')
>>> data = exareader(fnames)
>>> data.shape
(2733750,)   # this represents 1.5GB of simulation data
>>> data.dtype.names
('time', 'x', 'y', 'z', 'Q0', 'Q1', 'Q2', 'Q3', 'Q4')
>>> player = Play1D(data, "Q0 Q1 Q4".split(), "d sx e".split())
>>> plt.ion()
>>> player.createPlot()
>>> plt.title("My Fancy Simulation")
>>> # now resize plot, etc.
>>> plt.createMovie("mymovie.mpg")

After creating the movie, you can create more movies or invoke
"webmize mymovie.mpg" to quickly generate a mymovie.webm file in
the same folder to distribute over the web.
"""

from __future__ import print_function
import sys, argparse, logging

import numpy as np
import matplotlib.pyplot as plt

# TODO: cf. headless features at convergence analysis helpers.py
	#import matplotlib
	#matplotlib.use('Agg') # if headless
	#import matplotlib.pyplot as plt

from exahelpers import find_nearest, unzip, middle, first, last, cleandoc
from exahelpers import ExaFrontend
from exareader import ExaReader


# Quantities: name of conserved variable in ExaHypE -[mappedTo]-> human readable name
# TODO: Should be get from configfile or parameters.
default_quantities = ["Q0", "Q1", "Q4"]
default_quantities_readable = ["d", "sx", "e"]

logger = logging.getLogger("exaplayer")

class Player:
	"""
	The Player comes in two versions: A 2D plot or a 3D plot with color
	encoding. In any case, the output is a movie for the whole simulation
	time.
	"""
	def __init__(self, sol, quantities, names=None, zslice=None):
		self.sol = sol
		self.quantities = quantities
		self.names = names if names else quantities

		self.uni = lambda col: np.unique(sol[col])

		self.times = self.uni('time')
		self.x, self.y, self.z = map(self.uni, "xyz")
		self.spacedim = 3 if len(self.z) > 1 else 2

		self.zslice = zslice if zslice else middle(self.z)
	
	def something(self):
		raise NotImplementedError("Please Implement this method")

	@classmethod
	def add_group(cls, argparser):
		group = argparser.add_argument_group('plotting', description=cleandoc(cls))

		group.add_argument('-d', '--dimensions', choices=[1,2], type=int, default=2,
			           help="Whether to create movies of 1D or 2D plots")
		group.add_argument('-m', '--movie', metavar='somewhere.mpg', default=None,
			           help='Where to store a movie. If not given, no movie is stored')
		group.add_argument('-s', '--stay', action='store_true', default=False,
			           help='Whether to keep the window open after plotting, for interactivity.')
		group.add_argument('-I', '--interactive', action='store_true', default=True,
		                   help='Whether to switch on the pylab interactive flag')

	@staticmethod
	def main(args=None):
		"""
		Just call this as Player().main(sys.argv) or omit the args
		"""
		
		frontend = ExaFrontend(epilog=__doc__,
			program_description='ExaHyPE simulation data player')

		reader = ExaReader(files_as_main=False) # files_as_main=True does not work for me :(
		frontend.add_module(reader)
		frontend.add_module(Player)
		# TODO: slicer io_group could be used here to specify how to slice.
		# Actually this is now written at exaslicer.py.

		args = frontend.parse_args(args)

		logger.info("Welcome to vtkplayer, your friendly data plotter for ExaHyPE")
		#args = reader.remove_prog_in_inputfiles_list(frontend.parser, args)
		sol = reader.read_files_as_requested(args)
		registry = { 1: Play1D, 2: Play2D }
		PlayerClass = registry[args.dimensions]

		if args.interactive: plt.ion();
		logger.info("vtkplayer: Will do a %d-dimensial movie" % args.dimensions)

		player = PlayerClass(sol, default_quantities, default_quantities_readable)

		logger.info("vtkplayer: Having loaded %d timesteps in the range (%f,%f) in a %d-dimensional problem" % \
			(len(player.times), min(player.times), max(player.times), player.spacedim))

		player.createPlot()

		if args.movie:
			logger.info("Creating Movie...")
			player.createMovie(args.movie)

		if args.stay:
			logger.info("Calling Pyplots show, press CTRL+C to stop the program")
			plt.show(block=True)
		# if not stay:
		logger.debug("You could now start to interact with the code in IPython")

	def createMovie(self, output_filename, fps=15):
		"""
		Use this by hand to create movies after fixing the window size, setting
		titles and labes, etc. or use it from command line
		"""
		import matplotlib.animation as mani
		writers_i_like = ['ffmpeg', 'avconv']
		writers = [mani.writers[x] for x in writers_i_like if mani.writers.is_available(x)]
		if not writers: raise ValueError("Error: No writers in %s which I like" % mani.writers.avail)

		Writer = writers[0]
		#metadata = dict(title='Slice of %s' % quantities_readable, artist='ExaHyPE/vtkplayer')
		writer = Writer(fps=fps)#, metadata=metadata)
		logger.info("Starting movie output to file %s with writer %s" % (output_filename, writer))
		with writer.saving(plt.gcf(), output_filename, 100):
			for t in self.times:
				self.timeslider.set_val(t)
				plt.draw()
				writer.grab_frame()
		logger.info("Done. Movie at " + output_filename)


class Play1D(Player):
	def __init__(self, sol, quantities, names=None, zslice=None, yslice=None):
		"""
		Play1D requires same arguments as Player, but also
		a yslice. If not given, takes middle(y).
		"""
		Player.__init__(self, sol, quantities, names, zslice)
		self.yslice = yslice if yslice else middle(self.y)

	def slice1d(self, time, quantity):
		data1d = self.sol[
			(self.sol['time'] == time) &
			(self.sol['y'] == self.yslice) &
			(self.sol['z'] == self.zslice)
		]
		sx, ix = np.unique(data1d['x'], return_index=True)
		sy = data1d[quantity][ix]
		return (sx,sy)

	def createSubPlot(self, ax, quantity):
		sx, sy = self.slice1d(self.times[0], quantity=quantity)
		line, = ax.plot(sx, sy, "o-")
		title = ax.set_title("title")
		ax.set_xlim(0,1) # ExaHyPE coordinate system range
		return (line,title)

	def draw_update(self, slider_value=0):
		time = find_nearest(self.times, slider_value)
		logger.info("Plot for time="+str(time))
	
		for i, qu in enumerate(self.quantities):
			tt = "Slice of %s at t=%f, y=%.1f, z=%.1f" % (self.names[i], time, self.yslice, self.zslice )
			self.axes[i].set_title(tt)
			self.lines[i].set_data(self.slice1d(time, quantity=qu))
			plt.show()

	def createPlot(self):
		self.gs = plt.GridSpec(len(self.quantities),1)
		self.axes = map(plt.subplot, self.gs)
		self.lines, titles = unzip([ self.createSubPlot(a,q) for a,q in zip(self.axes, self.quantities)])
		plt.subplots_adjust(left=0.10, top=0.85)
		axcolor = 'lightgoldenrodyellow'
		self.ax1 = plt.axes([0.10, 0.90, 0.8, 0.03], axisbg=axcolor)

		self.timeslider = plt.Slider(self.ax1, "Time", min(self.times), max(self.times), min(self.times))
		self.timeslider.on_changed(self.draw_update)

		plt.draw()

class Play2D(Player):
	def __init__(self, sol, quantities, names=None, zslice=None):
		"""
		Play2D requires same arguments as Player.
		"""
		Player.__init__(self, sol, quantities, names, zslice)
		self.X, self.Y = np.meshgrid(self.x, self.y)

		# also convert the solutions to a continuous array
		# and store the column names for later access
		self.con = sol.view(sol.dtype[0])
		self.fields = { f: int(i) for i,f in enumerate(sol.dtype.names) }
		# fields is like {'x':0, 'y':1, 'z':2, 'time':3, 'Q0':4, ...}

		self.ix = last(np.unique(self.sol['x'], return_inverse=True))
		self.iy = last(np.unique(self.sol['y'], return_inverse=True))
		# another variant: Niftily store more columns. However this
		# copies the array, with several GB not working.
		#
		# from numpy.lib.recfunctions import append_fields
		# append_fields(self.sol, ['ix','iy'], [ix,iy]

		# a mesh which has the structure as 2D should have it
		self.X, self.Y = np.meshgrid(self.x, self.y)

	def conField(self, fieldname):
		return self.con[:, self.fields[fieldname]]
		
	def better_slice2d(self, time, quantity):
		assert time in self.times, "Time %f not a valid time" % time
		data2d = self.con[ (self.conField('time') == time) & (self.conField('z') == self.zslice) ]
		if not len(data2d):
			raise IndexError("No data for time %f" % time)
		# unique X and Y (like self.X, self.Y) and their positions in the
		# data2d array
		X, ix = np.unique(data2d[:,self.fields['x']], return_index=True)
		Y, iy = np.unique(data2d[:,self.fields['x']], return_index=True)
		
		mesh = zip(ix,iy)

	def slice2d(self, time, quantity):
		assert time in self.times, "Time %f not a valid time" % time
		data2d = self.sol[ (self.sol['time'] == time) & (self.sol['z'] == self.zslice) ]
		if not len(data2d):
			raise IndexError("No data for time %f" % time)
		Z = np.zeros_like(self.X)

		# due to issue #50, this will copy same values over and over.
		for i, col in enumerate(data2d):
			Z[ self.ix[i], self.iy[i] ] = col[quantity]

		return Z


	def draw_update(self, slider_value=0):
		time = find_nearest(self.times, slider_value)
		logger.info("2D Plot for time=%f" % time)
		self.ax1.set_title("2D plots at t=%f and z=%f" % (time, self.zslice))
	
		for i, qu in enumerate(self.quantities):
			self.images[i].set_data(np.transpose(self.slice2d(time, quantity=qu)))
		plt.show()

	def createSubPlot(self, ax, quantity, quantity_name):
		InitialData = self.slice2d(self.times[0], quantity=quantity)
		image = ax.imshow(InitialData, interpolation='none', aspect='auto', \
			origin='lower', extent=[min(self.x),max(self.x),min(self.y),max(self.y)])
		title = ax.set_title("Slice of %s (%s)" % (quantity_name, quantity))
		plt.colorbar(image, ax=ax)
		ax.set_ylabel("y"); ax.set_xlabel("x")
		return image

	def createPlot(self):
		self.gs = plt.GridSpec(1, len(self.quantities))
		self.axes = map(plt.subplot, self.gs)
		self.images = [self.createSubPlot(*o) for o in zip(self.axes, self.quantities, self.names)]
		plt.subplots_adjust(left=0.10, top=0.85)

		axcolor = 'lightgoldenrodyellow'
		self.ax1 = plt.axes([0.10, 0.90, 0.8, 0.03], axisbg=axcolor)
		self.timeslider = plt.Slider(self.ax1, "Time", min(self.times), max(self.times), min(self.times))
		self.timeslider.on_changed(self.draw_update)

		plt.draw()


if __name__ == "__main__":
	Player.main() # sys.argv


