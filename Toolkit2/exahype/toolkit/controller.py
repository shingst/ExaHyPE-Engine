#!/usr/bin/env python3

"""
This mimics the classical Frontend of the ExaHyPE toolkit.
"""

import os, sys, argparse, subprocess, datetime
from os.path import isdir, isfile
from pathlib import Path
from collections import OrderedDict

sys.path.append(os.path.join(os.path.dirname(__file__),"..","..")) # to allow import exahype... work
from exahype.toolkit import *
from exahype.specfiles import validate

from exahype.toolkit.makefile import *


class Controller():
    def header(self):
        info = {
            "gittag": subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip().decode('ascii'),
            "year": str(datetime.datetime.now().year)
        }
        return """
    ______           __  __      ____  ______    
   / ____/  ______ _/ / / /_  __/ __ \/ ____/    *************************
  / __/ | |/_/ __ `/ /_/ / / / / /_/ / __/        The ExaHyPE Toolkit v2
 / /____>  </ /_/ / __  / /_/ / ____/ /___           www.exahype.eu
/_____/_/|_|\__,_/_/ /_/\__, /_/   /_____/          Commit: %(gittag)s
                       /____/                    *************************

 This project has received funding from the European Union's Horizon 2020
 research and innovation programme under grant agreement No 671698 (ExaHyPE).

 Copyright (c) %(year)s, All rights reserved
 ExaHyPE is based  on the PDE framework Peano (www.peano-framework.org).
 Released under the BSD 3 Open Source License.
 
""" % info
        # end of header

    def wait_interactive(self, msg):
        if self.interactive:
            print("\n\n\n\n")
            print(msg + "... ok")
            input("<pres Enter>")
            print("\n\n\n\n")
            
    def info(self, msg):
        if self.verbose:
            print(msg)
            
    def load(self, specfile_name):
        """
        Given a specfile name, it assumes it to be JSON, reads it, validates it
        against the JSON-Schema and returns the native python data structure,
        enriched with default values from the schema.
        """
        specification = validate(specfile_name)
        return specification
            
    def __init__(self):
        parser = argparse.ArgumentParser(
            description="ExaHyPE's new python-based toolkit",
            epilog="See http://www.exahype.eu and the Guidebook for more help."
        )
        
        # some mandatory arguments from Toolkit1
        parser.add_argument("-c", "--clean-opt",
            help="Clean optimized kernels (only applicable if optimized kernels are used in Specfile)")
        g = parser.add_mutually_exclusive_group()
        g.add_argument("-i", "--interactive", action="store_true", default=False, help="Run interactively")
        g.add_argument("-n", "--not-interactive", action="store_true", default=True, help="Run non-interactive in non-stop-mode (default)")
        parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose (off by default, triggered on by --interactive)")

        parser.add_argument('specfile',
            type=argparse.FileType('r'),
            help="The specification file to work on (can be .exahype, .exahype2, .json)")
        
        args = parser.parse_args()
        
        self.interactive = args.interactive or not args.not_interactive
        self.verbose = args.verbose or self.interactive
        
        if self.verbose: # otherwise no need to call git, etc. 
            self.info(self.header())
        
        inpath = Path(args.specfile.name)
        if inpath.exists():
            self.info("Read input file %s." % inpath.resolve())
        else:    self.info("Read from stdin (%s)" % str(args.specfile))

        try:
            spec = self.load(args.specfile)
        except Exception as e:
            print("Could not properly read specfile")
            print(e)
            sys.exit(-3)
        
        try:
            d = DirectoryAndPathChecker(spec, self.verbose)
        except BadSpecificationFile as e:
            print("Some directories did not exist and/or could not be created.")
            print("Error message: ", e)
            print("Failure due to bad specificaiton file, cannot continue")
            sys.exit(-4)
        
        self.wait_interactive("validated and configured pathes")
        
        try:
            s = SolverGenerator(spec, args.specfile.name, self.verbose)
            s.generate_all_solvers()
        except BadSpecificationFile as e:
            print("Could not create applications solver classes.")
            print(e)
            sys.exit(-6)
        
        self.wait_interactive("generated application-specific solver classes")
        
        try:
            # kernel calls
            k = KernelCallsGenerator(spec, args.specfile.name, self.verbose)
            k.generate_solver_registration()
        except BadSpecificationFile as e:
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
        
        try:
            # kernel calls
            m = MakefileGenerator(spec, args.specfile.name, self.verbose)
            m.generate_makefile()
        except BadSpecificationFile as e:
            print("Could not create application-specific makefile")
            print(e)
            sys.exit(-10)
        
        self.wait_interactive("generated application-specific makefile")
