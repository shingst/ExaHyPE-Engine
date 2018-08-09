#!/usr/bin/env python3

import sys
import os
from controller import Controller

def main():
    """"Check python version, call the controller and run it"""
    # check version. Python 3.3 required
    requiredVersion = (3,3)
    currentVersion  = sys.version_info
    if(requiredVersion > currentVersion):
        sys.exit("Requires Python 3.3 or newer. Abort.")
    
    # create controller and run it (input parsed with argparse)
    control = Controller()


if __name__ == "__main__":
    # execute only if run as a script
    main()

