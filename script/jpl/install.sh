#! /usr/bin/bash

# exit when any command fails
set -e

if [[ $EUID -ne 0 ]]; then
   echo "This script must be run as root"
   exit 1
fi

if ! test -d include ; then
  echo "Failed to find the \'include\' directory (under pwd)"
  exit 1
fi

if ! test -d lib ; then
  echo "Failed to find the \'lib\' directory (under pwd)"
  exit 1
fi

pwdd=$(pwd)
INC=${pwdd}/include
LIB=${pwdd}/lib

INC_DIR=/usr/local/include/cppspice
if ! test -d ${INC_DIR} ; then mkdir ${INC_DIR} ; fi

cd include
echo "Linking/Installing header files to ${INC_DIR}"
for hf in *.h ; do
  if test -f ${INC_DIR}/${hf} ; then rm ${INC_DIR}/${hf} ; fi
  ln -s ${INC}/${hf} ${INC_DIR}/${hf} 
done
cd ../


LIB_DIR=/usr/local/lib

cd lib
echo "Linking/Installing library files to ${LIB_DIR}"
for lf in *.a ; do
  if test -f ${LIB_DIR}/${lf} ; then rm ${LIB_DIR}/${lf} ; fi
  ln -s ${LIB}/${lf} ${LIB_DIR}/lib${lf} 
done
cd ../
