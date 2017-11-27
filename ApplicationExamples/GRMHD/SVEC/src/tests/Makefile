SHELL=/bin/bash

TEST=TEST_NEW_PDE_AUTONOMOUSLY

# path to GRMHD base directory
GRMHD=../../..

# for tarch etc
PEANO_INC=../../../../../Peano
EXAHYPE_INC=../../../../../ExaHyPE

# GRMHD Fortran formulation
#GRMHD_Finc=$(GRMHD)/GRMHD
FortranObjects=$(GRMHD)/GRMHD/Fortran
FortranInitialData=$(GRMHD)/GRMHD/InitialData

# GRMHD cpp formulation
GRMHD_Cinc=$(GRMHD)/GRMHD_cpp

# # WAS Dim3
DIMSTRING=Dim2

all:
	g++ -c -g3 --std=c++11 -D$(DIMSTRING) -D$(TEST) -I$(GRMHD_Cinc) -I../ -I$(PEANO_INC) -I$(EXAHYPE_INC) -Wall \
		indi-test.cc $(GRMHD_Cinc)/InitialData/AnalyticalID.cpp \
		../PDE/Cons2Prim.cpp ../PDE/Prim2Cons.cpp ../PDE/PDE-GRMHD.cpp

#	We link the original GRMHD (Fortran) files here. Compile them by running:
#
#	cd ../../../GRMHD/
#	[ -e Makefile ] || exa toolkit GRMHD_ADERDG.exahype
#	for x in {InitialData,Fortran}/*.{f90,cpp}; do sh -c "make $(pwd)/${x%.*}.o"; done;

		
	g++ -g3 -oinditest.out \
		$(FortranInitialData)/InitialDataFort.o \
		$(FortranObjects)/{abort,C2P-GRMHD,C2PRoutines,Metric,Parameters,PDE}.o \
		*.o -lgfortran

clean:
	rm -f *.o *.out
