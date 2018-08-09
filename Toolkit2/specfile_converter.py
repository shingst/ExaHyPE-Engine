#!/usr/bin/env python3
# run e.g.: ./specfile_converter.py ../Demonstrators/EulerADERDG/EulerADERDG.exahype
import sys
import configparser
import re

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
solver=0
plotter=0
spec_file_1_cleaned = ""
reads_multiline_comment = 0
reads_singline_comment  = False
for line in spec_file_1:
	reads_singline_comment = line.strip().startswith("//") or line.strip().startswith("#")
	reads_multiline_comment += 1 if line.strip().startswith("/*") else 0
	if reads_multiline_comment==0 and not reads_singline_comment:
		m_project = re.match(r"^\s*exahype-project\s+(\w+)",line)
		m_group   = re.match(r"^\s*(computational-domain|shared-memory|distributed-memory|global-optimisation)",line)
		m_solver  = re.match(r"^\s*solver\s+([^\s]+)\s+(\w+)",line)
		m_plotter = re.match(r"^\s*plot\s+([^\s]+)\s+(\w+)",line)
		if m_project:
			spec_file_1_cleaned += "[exahype-project]\n"
			spec_file_1_cleaned += "name="+m_project.group(1)+"\n"
		elif m_group:
			spec_file_1_cleaned += "[%s]" % m_group.group(1)+"\n"
		elif m_solver:
			spec_file_1_cleaned += "[solver%d]" % (solver)  +"\n"
			spec_file_1_cleaned += "solver_type="+m_solver.group(1)+"\n"
			spec_file_1_cleaned += "name="+m_solver.group(2)+"\n"
			solver+=1
		elif re.match(r"end\s*solver",line):
			plotter=0
		elif m_plotter:
			spec_file_1_cleaned += "[solver%dplotter%d]\n" % (solver,plotter)
			spec_file_1_cleaned += "plotter_type="+m_plotter.group(1)+"\n"
			spec_file_1_cleaned += "name="+m_plotter.group(2)+"\n"	
			plotter+=1
		elif line.strip().startswith("end"):
			pass
		else:
			spec_file_1_cleaned += line.strip().replace("const","")+"\n"
	reads_multiline_comment -= 1 if line.strip().startswith("*/") else 0 

print(spec_file_1_cleaned)

config = configparser.ConfigParser()
config.read_string(spec_file_1_cleaned)

print(config.sections())

# now we need the mapping
