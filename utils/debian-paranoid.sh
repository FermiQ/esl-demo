#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.0.1
mkdir build
cd build
export FC=gfortran
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -Wall -fcheck=all -finit-integer=snan -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow " cmake ../
make -j VERBOSE=1
