#!/bin/bash

# Would be better if there was a way to know which one is the latest snapshot...

# Will remove local Peano installation

# arg: download the n-th snapshoot found, if no argument then the first
if [[ $1 =  "" ]];
then
  n=1
else
  n=$1
fi

# Make sure the script can be called from anywhere.
cd "$(dirname "$0")"
cd ../Submodules

# Download page
URL=http://www.peano-framework.org/download/peano/

# Regex to match snapshoot
re="\"(Peano-.{7}\.tar\.gz)\""

match[0]=""
global_rematch() { 
    local s=$1 regex=$2 i=1
    while [[ $s =~ $regex ]]; do 
        match[$i]="${BASH_REMATCH[1]}"
        s=${s#*"${BASH_REMATCH[1]}"}
        i=$i+1
    done
}


# Get HTML
content=$(curl -L $URL)
# Find snapshoot name
global_rematch "$content" "$re" 

if [[ $content =~ $re ]]; 
then
	FILE_NAME=${match[$n]}
  if [[ $FILE_NAME =  "" ]];
  then
    echo "Couldn't find the source, check $URL"
  else
    echo "Snapshoot found: $FILE_NAME"
    DOWNLOAD_URL=$URL''$FILE_NAME
  	echo "Download snapshoot"
  	wget $DOWNLOAD_URL
  	echo "Remove existing installation"
  	rm -rf Peano
  	echo "Extract snapshoot"
  	tar -xf $FILE_NAME src/
    mkdir Peano
  	mv src Peano
  	echo "Cleanup"
  	rm -rf *tar.gz
  fi
else
	echo "Couldn't find the source, check $URL"
fi
