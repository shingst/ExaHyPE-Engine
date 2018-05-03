SHELL=/bin/bash

# for tarch etc
EXAROOT=../../..
PEANO_INC=$(EXAROOT)/Peano
EXAHYPE_INC=$(EXAROOT)/ExaHyPE
# Dim2 or Dim3
DIMSTRING=Dim2

# This is a naive way to link to the existing object files of an ExaHyPE installation.
# It requires a compiled application with the same dimensionality.
# You may have linking errors from old ("dangling") object files where no more cpp file
# exists to. Delete these object files as a make clean in the application does not
# catch them.
OFILES =$(shell find $(PEANO_INC)/ -iname \*.o)
OFILES+=$(shell find $(EXAHYPE_INC)/ -iname \*.o)

# remove some files again
OFILES := $(filter-out $(EXAHYPE_INC)/exahype/main.o,$(OFILES))
OFILES := $(filter-out $(EXAHYPE_INC)/exahype/Version.o,$(OFILES))


parsertest: libexahype.a
	g++ -c -g3 --std=c++11 -D$(DIMSTRING) -I$(PEANO_INC) -I$(EXAHYPE_INC) -Wall \
		test-parser.cc
	g++ -c -g3 --std=c++11 mexa.cpp

	g++ -O0 -otest-parser.out test-parser.o mexa.o libexahype.a

libexahype.a:
	# All object files as an archive. This speeds up linking enormously.
	ar cr libexahype.a $(OFILES)
	
mexatest:
	g++ -std=c++11 mexa.cpp -Wall -DMEXA_TEST_MAIN -osimple-mexa.out

clean:
	rm -f *.a *.o *.out
