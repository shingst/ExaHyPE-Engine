#!/usr/bin/env python3
import sys
import configparser

def parse_argument(i):
	if i<len(sys.argv):
		return sys.argv[i]
	else:
		return None

spec_file_1 = ""
try:
	with open(parse_argument(1), "r") as input_file:
		spec_file_1 = input_file.readlines()
except Exception:
	raise

## remove comments
#print(spec_file_1
solver=-1
plotter=-1
spec_file_1_cleaned = ""
reads_multiline_comment = 0
reads_singline_comment  = False
for line in spec_file_1:
	reads_singline_comment = line.strip().startswith("//") or line.strip().startswith("#")
	reads_multiline_comment += 1 if line.strip().startswith("/*") else 0
	if reads_multiline_comment==0 and not reads_singline_comment:
		if line.strip().startswith("solver"):
			parts = line.split(" ")
			solver+=1
		if line.strip().startswith("end solver"):
			solver+=0
		if line.strip().startswith("plotter"):
			solver+=1
		if line.strip().sswith("plotter"):
			solver+=-1

	reads_multiline_comment -= 1 if line.strip().startswith("*/") else 0 

print(spec_file_1_cleaned)

config = configparser.ConfigParser()
