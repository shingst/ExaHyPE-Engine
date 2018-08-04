#!/usr/bin/env python3

"""
This mimics the classical Frontend of the ExaHyPE toolkit.
"""

import os, sys, argparse
from os.path import isdir, isfile
from pathlib import Path
from collections import OrderedDict

from . import BadSpecificationFile

from directories import DirectoryAndPathChecker
from solver import SolverClassGenerator
from plotter import PlotterGenerator
from kernelcalls import KernelCallsGenerator

class ClassicFrontend():

	def wait_interactive(self, msg):
		if self.interactive:
			print("\n\n\n\n")
			print(msg + "... ok")
			input("<pres Enter>")
			print("\n\n\n\n")
			
	def info(self, msg):
		if self.verbose:
			print(msg)

	def main(self):
		parser = argparse.ArgumentParser(
			description="ExaHyPE's new python-based toolkit",
			epilog="See http://www.exahype.eu and the Guidebook for more help."
		)
		
		# some mandatory arguments from Toolkit1
		parser.add_argument("-c", "--clean-opt",
			help="Clean optimized kernels (only applicable if optimized kernels are used in Specfile)")
		parser.add_argument("-i", "--interactive", help="Run interactively") # TODO: entangle with not-interactive by ArgumentParser
		parser.add_argument("-n", "--not-interactive", help="Run non-interactive in non-stop-mode (default)")
		parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose (off by default, triggered on by --interactive)")

		parser.add_argument('specfile',
			type=argparse.FileType('r'),
			help="The specification file to work on (can be .exahype, .exahype2, .json)")

		args = parser.parse_args()

		self.interactive = args.interactive or not args.not_interactive
		self.verbose = args.verbose or interactive
		
		self.info("Read input file %s." % Path(args.specfile).resolve())

		try:
			spec = parse_the_specfile(args.specfile) # TODO: Implement
		except e:
			print("Could not properly read specfile")
			print(e)
			sys.exit(-3)
		
		try:
			d = DirectoryAndPathChecker(spec, verbose)
		except BadSpecificationFile e:
			print("Some directories did not exist and/or could not be created.")
			print(e)
			sys.exit(-4)
		
		self.wait_interactive("validated and configured pathes")
		
		try:
			SolverGenerator(spec, verbose)
		except BadSpecificationFile e:
			print("Could not create applications solver classes.")
			print(e)
			sys.exit(-6)
	
		self.wait_interactive("generated application-specific solver classes")


		# This is no more possible:
		# We need to generate the plotters straight when looping over the Solvers,
		# since otherwise we lack information to which solvers they belong to
		#
		#try:
		#	PlotterGenerator(spec, verbose)
		#except BadSpecificationFile e:
		#	print("Could not create application's plotter classes")
		#	print(e)
		#	sys.exit(-8)
		#	
		#self.wait_interactive("generate application-specific plotter classes")
		
		try:
			# kernel calls
			KernelCallsGenerator(spec, verbose)
		except BadSpecificationFile e:
			print("Could not create ExaHyPE's kernel calls")
			print(e)
			sys.exit(-10)
			
		self.wait_interactive("generated computational kernel calls")
		
		# blah
		#CreateReadme(spec,verbose)
		#
		#self.wait_interactive("generated README")
		#
		# makefiles, etc.
		#setupBuildEnvironment(spec, verbose)
		#
		#self.wait_interactive("setup build environment")
