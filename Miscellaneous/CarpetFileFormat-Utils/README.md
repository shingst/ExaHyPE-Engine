# Utilities to handle CarpetHDF5 files

This directory contains standalone executables which allow conversion and
modification of CarpetHDF5 files. They were copied together from various
sources from Cactus thorns.

They were collected by Sven K at 2017-06-06.

## Example invocations

```
hdf5_merge output.time*.h5 output-single.h5
hdf5_slicer --out2d-xyplane-z 0 data.file*.h5 data.z=0.h5
hdf5_double_to_single data-double64bit.h5 data-float32bit.h5
```

## Sources

The files

```
hdf5_recombiner.cc
hdf5_slicer.cc
hdf5toascii_slicer.cc
hdf5tobinary_slicer.cc
```

in this directory stem from the Carpet/CarpetIOHDF5/src/utils/ directory, public
available as Open Source at https://bitbucket.org/eschnett/carpet

The files

```
hdf5_double_to_single.c
hdf5_extract.c
hdf5_merge.c
```

stem from the Cactus HDF5 thorn (`arrangements/ExternalLibraries/HDF5/src/util`),
public available as Open Source at https://svn.cactuscode.org/projects/ExternalLibraries/HDF5/

