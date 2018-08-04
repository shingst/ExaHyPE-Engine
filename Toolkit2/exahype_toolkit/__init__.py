#!/usr/bin/env python3

from . import frontend

run = frontend.main

class BadSpecificationFile(Exception):
	pass

if __name__ == '__main__':
	frontend.main()
