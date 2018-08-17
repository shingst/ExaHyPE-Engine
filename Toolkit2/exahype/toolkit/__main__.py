#!/usr/bin/env python3

import sys
import os

def main():
    """"Call the controller and run it"""
    # append directory above to sys.path to be able to load the module
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")))
    # import as module
    from toolkit import Controller
    control = Controller()
    control.run()


if __name__ == "__main__":
    # execute only if run as a script
    main()

