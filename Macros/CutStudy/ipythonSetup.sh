#!/bin/bash

source /afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/setup.sh #gcc compiler

pwd=$PWD
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.24/x86_64-slc6-gcc47-opt/root
source bin/thisroot.sh
cd ${pwd}

echo "gcc set to:"
which gcc

echo " "
echo "root version set to:"
echo $ROOTSYS
