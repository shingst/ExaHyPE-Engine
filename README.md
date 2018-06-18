# This is the ExaHyPE project #

This is the maing repository of the ExaHyPE project. It contains the following directories:

* ApplicationExamples: A number of exemplaric application such as Euler, MHD, Seismic, etc.
* Benchmarks: An ExaHyPE Benchmark suite for TBB and MPI (maintained by Durham)
* CodeGenerator: The optimized kernels generators (maintained by TUM)
* Demonstrators: Lightweight applications for demonstration purpose (maintained by Durham)
* Miscellaneous: Various codes and tools
* Peano: Mounting point for the Peano installation
* Toolkit: The ExaHyPE development toolkit
* UserCodes: Prototypes of ExaHyPE written in Matlab and Fortran and contributed by Trento

## Mini installation guide ##

Copy and paste these commands to start with a working ExaHyPE application and compile the demo application _EulerFlow_:

    git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
    ./ExaHyPE-Engine/Peano/checkout-update-peano.sh
    ./ExaHyPE-Engine/Toolkit/build.sh

Now you are ready to follow compile and run an ExaHyPE application [according to the guidebook](http://www5.in.tum.de/exahype/guidebook.pdf):

    cd ExaHyPE-Engine/
    java -jar Toolkit/dist/ExaHyPE.jar ApplicationExamples/EulerFlow.exahype
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

### CodeGenerator dependencies (optional) ###

The optimised kernels are optional and can be replaced by the generic ones, thus the CodeGenerator's dependencies are also optional.

Python 3 is required to run the CodeGenerator.

Additionally the CodeGenerator requires:

* libxsmm (https://github.com/hfp/libxsmm) to generate the advanced matrix multiplication code using libxsmm's generator
* Jinja2 (https://github.com/pallets/jinja.git) a python3 template engine to generate the optimised kernel
* MarkupSafe (https://github.com/pallets/markupsafe.git), a dependency from Jinja2.

A script is provided to import all the dependencies locally.
In contrast to the CodeGenerator, the build of libxsmm 
requires Python 2. It fails with Python 3.
So, make sure that "python" does refer
to a Python 2 binary.

Quick installation:

    ./ExaHyPE-Engine/CodeGenerator/importDependenciesLocally.sh 

If Jinja2 is provided by your python installation you can safelly edit the value of the import script configuration parameter ``JINJA2_ALREADY_AVAILABLE`` to remove its and MarkupSafe local imports.


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



