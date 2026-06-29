#!/bin/bash

mkdir -p  $PREFIX/share/crpropa &> /dev/null  # eliminate possible warning outputs
cd $PREFIX/share/crpropa
CRPROPA_DATAFILE_VER=$(cat CMakeLists.txt | grep --extended-regexp "set\( *CRPROPA_DATAFILE_VER +.+\)" | sed -E "s/set\( *CRPROPA_DATAFILE_VER +\"//g" | sed -E "s/\"\)//g")
USERPWD=3juW9sntQX2IWBS

# download data
curl -u $USERPWD:$USERPWD -O https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
curl -u $USERPWD:$USERPWD -O https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz

# check download hash, command depends on OS:
OS=$(uname)
if [ $OS = "Linux" ]
then
	echo $(cat data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM) | md5sum --check --status data-${CRPROPA_DATAFILE_VER}.tar.gz
elif [ $OS = "Darwin" ]
then
	if [ $(cat data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM | cut -d ' ' -f 1) = $(md5 data-${CRPROPA_DATAFILE_VER}.tar.gz | cut -d ' ' -f 4) ]
	then
		:
	else
		exit 1
	fi
fi

tar xzf data-${CRPROPA_DATAFILE_VER}.tar.gz
mv data/* $PREFIX/share/crpropa/
# cleanup
rm -rf data data-${CRPROPA_DATAFILE_VER}.tar.gz data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
