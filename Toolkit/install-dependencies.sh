#!/bin/bash
#
# This is an interactive script to guide throught the dependencies installation
# for the Python3 packages needed to run the ExaHyPE toolkit. It is safe to
# call as it will ask before doing anything.
#

# TODO SVEN: your install script won't work, it need to install to Submodules or change the paths in the Configuration files of toolkit, specfiles and codegenerator
echo "TODO SVEN: your install script won't work, it need to install to Submodules or change the paths in the Configuration files of toolkit, specfiles and codegenerator"
exit -1


Toolkit2="$(dirname $0)"
has() { type $@ &>/dev/null; } # a way to check if command is available

GITHUB=git://github.com

# Uncomment if you use port forwarding a la
# ssh -R 19489:github.com:9418 <Your Login>@<SuperMUC Login Node>
#GITHUB=git://localhost:19489

modules=(attr pyresistent markupsafe jinja2 jsonschema)

repos=(
  $GITHUB/python-attrs/attrs.git\
  $GITHUB/tobgu/pyrsistent.git\
  $GITHUB/pallets/markupsafe.git\
  $GITHUB/pallets/jinja.git\
  $GITHUB/Julian/jsonschema.git\
)

echo "$0: In order to install the ExaHyPE toolkit dependencies, you can"

if has pip3; then
  echo "$0: - Use the Python package manager 'pip' to install the packages:"
  echo "$0:"
  echo "$0:   sudo pip3 install ${modules[@]}   # system wide"
  echo "$0:   pip3 -U install ${modules[0]}     # only for user $(whoami)"
  echo "$0:"
fi

if has git; then
  echo "$0: - Download the codes from their git repositories, via $GITHUB,"
  echo "$0:   to a local subfolder (guided by this script)."
  echo "$0:"
  read -p "$0:   Do you want to install dependencies in local subfolder (y/n)?" yn
  case $yn in
    [Yy]* ) 
        dependencies=$Toolkit2/dependencies
        mkdir -vp $dependencies
        rm -rvf $dependencies/*
      
        for i in ${!modules[*]}; do
          module=${modules[$i]}
          repo=${repos[$i]}
                            
          # symlink dependencies into top level dir 
          ( cd $dependencies && git clone --verbose $repo ) # returns to work dir afterwards
          
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

        echo "$0: Finished successfully installing all dependencies."
      ;;
    #[Nn]* )
    * )
        echo "$0: Answer considered as no."
        exit -1
      ;;
  esac
else
  echo "$0: - NOT download the codes from their git repositories since you don't have the git repository on your path."
fi
