
SHELL=bash

CXX = h5c++
EXE = binns-dumper
TARGETS = hdf5_recombiner hdf5_slicer hdf5toascii_slicer hdf5tobinary_slicer hdf5_double_to_single hdf5_merge hdf5_extract

all: $(TARGETS)

$(TARGETS):
	$(CXX) -o $@ -O2 $@.cc

clean:
	rm -f $(TARGETS) *.o
