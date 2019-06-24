# This is the ExaHyPE project #

This is the main repository of the ExaHyPE project*. It contains the following directories:

* ApplicationExamples: A number of exemplaric application such as Euler, MHD, Seismic, etc.
* Benchmarks: An ExaHyPE Benchmark suite for TBB and MPI (maintained by Durham)
* CodeGenerator: The optimized kernels generators (maintained by TUM)
* Demonstrators: Lightweight applications for demonstration purpose (maintained by Durham)
* Miscellaneous: Various codes and tools
* Peano: Mounting point for the Peano installation
* Submodules: Installation directory for the git submodules (ExaHyPE's dependencies)
* Toolkit: The ExaHyPE development toolkit
* UserCodes: Prototypes of ExaHyPE written in Matlab and Fortran and contributed by Trento

## Mini installation guide ##

### Dependencies ###

Python3 (3.3 at least) is required for the Toolkit and the CodeGenerator (compile time).

Other dependencies are open source and registered as git submodules, a working C++ compiler is required to build libxsmm. 
A script is provided to fetch them, and if necessary build them.

### Instalation and demo ###

Copy and paste these commands to start with a working ExaHyPE application and compile the demo application _EulerFlow_:

    git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
    ./ExaHyPE-Engine/Submodules/updateSubmodules.sh

Now you are ready to follow compile and run an ExaHyPE application [according to the guidebook](http://www5.in.tum.de/exahype/guidebook.pdf):

    cd ExaHyPE-Engine/
    ./Toolkit/toolkit.sh ApplicationExamples/EulerFlow.exahype
    cd ApplicationExamples/EulerFlow/
    
    export COMPILER=gnu
    export TBB_INC=/usr/include/tbb
    export TBB_LIB=/usr/lib/tbb
    make
    ./ExaHyPE-Euler ../EulerFlow.exahype

Look into `RUN.sh` for an alternative, more elaborate way how to setp the installation and it's dependencies. You also might want to use the `exa` tool which makes it super simple to start and use an installation from the scratch: After downloading, just execute

    exa=Miscellaneous/BuildScripts/exa.sh 
    $exa bootstrap
    $exa compile-run EulerFlow
    

## General remarks ##

* Run tests before you commit
* Document your code with doxygen
* Disable auto-formatting of your IDE or follow the google code-style => `Code/Miscellaneous/.clang-format`. For Eclipse users, there is also https://github.com/wangzw/CppStyle.
* Do not run autoformatters on the DaStGen definition files (`*.def`). This will screw up the `Packed-type: ..` and `Constant: ..` lines.


## Commit guidelines ##

Please, don't commit the following:
    
* Binary files (`*.o, executables, ... `) excluding those necessary for the documentation 
* Output files (`*.vtk, logs, ... `)

Please write good commit messages that document how you changed ExaHyPE.


## Regenerate ExaHyPE's kernel gluecode ##
 
```
java -jar ~/workspace/peano/pdt/pdt.jar --generate-gluecode exahype/exahype.specification exahype ~/workspace/peano/pdt/usrtemplates:../Peano/multiscalelinkedcell
```


## Build a new release ##

Recommendation: create a release separate from your working copy, e.g. in a ExaHyPE release folder.

1) Clone/pull the ExaHyPE-Engine repository

2) Run `Peano/checkout-update-peano.sh`

3) Clone/pull the ExaHyPE-Documentation repository

4) Run `make release` in ExaHyPE-Documentation

5) Run `upload.sh` in ExaHyPE-Documentation
    - This uploads a copy of the guidebook to http://dev.exahype.eu/guidebook.pdf

6) Run `Miscellaneous/create-release.sh` and pass the application folder(s) you would like to add to the release

7) This will create a tar.gz archive containing the ExaHyPE+Peano source and your application(s)


* This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 671698.
