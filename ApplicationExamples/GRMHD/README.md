# The GRMHD equations in ExaHyPE

This directory holds five subdirectories:

 - doc: Latex documentation about the PDEs. Compile and read them.
 - GRMHD: The original F90 formulation of the GRMHD PDEs. They used to work
   in principle but had some bugs with C++ linking in the C2P.
 - GRMHD_f90: Really old variant which used the Fortran kernels.
 - SVEC: A tensor code for the C++ formulation of the GRMHD equations
 - GRMHD_cpp: The novel (end 2017) PDE written in pure C++

## How to compile the Hybrid applications

If you look into the respective application directories, you will notice the solvers:

```
GRMHDSolver_ADERDG.cpp
GRMHDSolver_ADERDG.h
GRMHDSolver_FV.cpp
GRMHDSolver_FV.h
```

In order to compile these files, you need to run the toolkit on both the FV and ADERDG
versions of that application. That is, you should do three runs, for instance to
compile the GRMHD_cpp in 3D:

  1. Toolkit GRMHD_cpp_FV-3D.exahype
  2. Toolkit GRMHD_cpp_ADERDG-3D.exahype
  3. Toolkit *the actual specfile you want to use*
  4. Compile

If you don't run the toolkit on both variants (FV,ADERDG), the files will not compile
because they complain about missing classes which are generated by the toolkit.

-- SK, 2017-10-24