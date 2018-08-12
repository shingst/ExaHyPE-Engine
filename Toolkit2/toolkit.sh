#!/bin/bash
#
# This is a neat frontend to the ExaHyPE toolkit NG.
# It checks whether all Python is available and then calls it.
# 

Toolkit2="$(dirname $0)"
has() { type $@ &>/dev/null; } # a way to check if command is available


# find suitable python on the path
if has python3; then PYTHON3="python3";
elif python --version | grep -qi "python 3"; then PYTHON3="python"
else echo "$0: Python3 required for running the ExaHyPE toolkit" >&2; exit -1; fi

# check that all required modules are there.
# Could probably postpone that because it is slow to call python so many times.


GITHUB=git://github.com

# Uncomment if you use port forwarding a la
# ssh -R 19489:github.com:9418 <Your Login>@<SuperMUC Login Node>
#GITHUB=git://localhost:19489

modules=(\
  attr\
  pyrsistent\
  markupsafe\
  jinja2\
  jsonschema\
  )

repos=(\
  $GITHUB/python-attrs/attrs.git\
  $GITHUB/tobgu/pyrsistent.git\
  $GITHUB/pallets/markupsafe.git\
  $GITHUB/pallets/jinja.git\
  $GITHUB/Julian/jsonschema.git\
  )


errors=false
for i in ${!modules[*]}; do
  module=${modules[$i]}
  if ! $PYTHON3 -c "import sys; sys.path.append(\"$Toolkit2\"); import $module" 2>&1 >/dev/null; then
    echo "$0: Required python3 module '$module', not available." >&2
    errors=true
  fi
done

#exec $PYTHON3 $Toolkit2/exahype/toolkit/frontend.py $@
if [ "$errors" = "false" ]; then
  PYTHONPATH="$Toolkit2" $PYTHON3 $Toolkit2/exahype/toolkit $@   # RUN PROGRAM
else
  echo "$0: There are missing dependencies." >&2
  echo "$0:"
  echo "$0: Either run as root: " >&2
  echo "$0:"
  echo "$0: pip3 install ${modules[@]}" >&2
  echo "$0:"
  echo "$0: Or install dependencies in local subfolder." >&2
  read -p "$0: Do you want to install dependencies in local subfolder (y/n)?" yn
  case $yn in
    [Yy]* ) 
        dependencies=$Toolkit2/dependencies
      
        if [ ! -d $dependencies ]; then
        mkdir $dependencies
      else
        rm -rf $dependencies/*
      fi
      
      for i in ${!modules[*]}; do
        module=${modules[$i]}
        repo=${repos[$i]}
                          
         # symlink dependencies into top level dir 
        ( cd $dependencies && git clone $repo ) # returns to work dir afterwards
        
        if [ "$module" == "attr" ]; then
          (cd $Toolkit2 && ln -sf dependencies/attrs/src/attr ./)
        elif [ "$module" == "jinja2" ]; then
          (cd $Toolkit2 && ln -sf dependencies/jinja/jinja2 ./)
        else
          (cd $Toolkit2 && ln -sf dependencies/$module/$module ./)
        fi
        # remove last two lines of jsonschema __init__.py file as module might not be registered when installed locally 
        # (head cannot use same file for input and output)
        if [ "$module" == "jsonschema" ]; then
          head -n -2 $Toolkit2/jsonschema/__init__.py > $Toolkit2/jsonschema/__init__.py.tmp  
          mv $Toolkit2/jsonschema/__init__.py.tmp $Toolkit2/jsonschema/__init__.py
        fi
      done
      ;;
    [Nn]* ) errors=true;;
    * ) errors=true; echo "Answer considered as 'n'.";;
  esac
fi
