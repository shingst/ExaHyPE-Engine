# This is the ExaHyPE project #

## Mini installation guide ##

Copy and paste these commands to start with a working ExaHyPE application:

    git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
    cd ExaHyPE-Engine/Code
    ./checkout-update-peano.sh
    cd ../Toolkit/
    ./build.sh
    cd ../
    ln -s Toolkit/dist/ExaHyPE.jar

Now you are ready to follow compile and run an ExaHyPE application [according to the guidebook](http://www5.in.tum.de/exahype/guidebook.pdf):

    java -jar ExaHyPE.jar ApplicationExamples/EulerFlow.exahype
    export COMPILER=gnu
    export TBB_INC=/usr/include/tbb
    export TBB_LIB=/usr/lib/tbb
    cd ApplicationExamples/EulerFlow && make -j
    ./ExaHyPE-Euler ../EulerFlow.exahype
    
### Libxsmm (optional) ###

Libxsmm's code generator is required by the Toolkit when using application taylored optimised kernels to generate the advanced matrix multiplication code. 
The optimised kernels are optional and can be replaced by the generic ones, thus libxsmm is also optional.

Libxsmm's sources can be found on github: https://github.com/hfp/libxsmm

Quick installation:

    cd ExaHyPE-Engine/
    git clone https://github.com/hfp/libxsmm.git
    cd libxsmm/
    make generator
    cd ../Code/
    ln -s ../libxsmm/ Libxsmm

When using the optimised kernels, the path to libxsmm has to be specified in the application's specification file.

## General remarks ##

* Run tests before you commit
* Document your code with doxygen
* Disable auto-formatting of your IDE or follow the google code-style => `Code/Miscellaneous/.clang-format`. For Eclipse users, there is also https://github.com/wangzw/CppStyle.
* Do not run autoformatters on the DaStGen definition files (`*.def`). This will screw up the `Packed-type: ..` and `Constant: ..` lines.


## Commit guidelines ##

Please, don't commit the following:
    
* Binary files (`*.o, executables, ... `) excluding those necessary for the documentation 
* Output files (`*.vtk, logs, ... `)

### Commit message template 
We should try to write good commit messages that document
every change we introduced into the code.

What do you think about the following template?  
(Source: http://www.git-scm.com/book/en/v2/Distributed-Git-Contributing-to-a-Project)

Template:  

    Short (50 chars or less) summary of changes

    More detailed explanatory text, if necessary.  Wrap it to
    about 72 characters or so.  In some contexts, the first
    line is treated as the subject of an email and the rest of
    the text as the body.  The blank line separating the
    summary from the body is critical (unless you omit the body
    entirely); tools like rebase can get confused if you run
    the two together.
    
    Further paragraphs come after blank lines.
    
    - Bullet points are okay, too
    
    - Typically a hyphen or asterisk is used for the bullet,
      preceded by a single space, with blank lines in
      between, but conventions vary here"

