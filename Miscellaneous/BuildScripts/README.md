# The exa build managament tools

This directory contains a number of scripts which are helpful to simplify the
build process of ExaHyPE. In the order of relevance, this is:

## The exa toolkit

Call ./exa.sh to uncover the swiss army knife of ExaHyPE build managament.

You should install it to your $PATH so you can invoke "exa" from any
directory. To do so, use ./install-exa.sh.

If you have to manage multiple installations on your account, use
install-relocatable-exa.sh

## Compile.sh

The ultimate compilation script is ./compile.sh. It bundles the environment
managament, toolkit invocation and make invocation including logfile generation.
It also manages to upload the build log to a pastebin in order to share it
with your collegues. It is solely controled via environment variables.

This script is the central workhorse in the whole exa build managament tools.

## Out of order compilation for ExaHyPE

If you need to run a parameter study with compile-time parameters in ExaHyPE
or you want to quickly compile many variants or applications in parallel, you
will have trouble with ExaHyPE's build system because it works *in-place*,
that means next to every .cpp file the compiled .o file is stored. What you
need is an *out-of-order* build system which seperates the code from the
compiled intermediate object files.

The setup-out-of-tree.sh script does this for you. Inspect the script for
further documentation.

