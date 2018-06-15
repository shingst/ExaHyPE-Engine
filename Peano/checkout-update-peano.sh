#!/bin/bash

# Make sure the script can be called from anywhere.
cd "$(dirname "$0")"

# Download page
URL=http://csdemo.ddns.net/download/peano/

# Regex to match snapshoot
re="\"(Peano-.{7}\.tar\.gz)\""

# Get HTML
content=$(curl -L $URL)
# Find snapshoot name
if [[ $content =~ $re ]]; 
then
	FILE_NAME=${BASH_REMATCH[1]} 
	echo "Snapshoot found"
	DOWNLOAD_URL=$URL''$FILE_NAME
	echo "Download snapshoot"
	wget $DOWNLOAD_URL
	echo "Remove existing installation"
	rm -rf doxygen-html peano tarch toolboxes
	echo "Extract snapshoot"
	tar -xf $FILE_NAME src/
	mv src/* .
	echo "Cleanup"
	rm -rf *tar.gz
	rm -rf src
else
	echo "Couldn't find the source, check $URL"
fi