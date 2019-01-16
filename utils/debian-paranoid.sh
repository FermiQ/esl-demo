#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.3.0
mkdir build-paranoid
cd build-paranoid
export FC=gfortran
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -Wall -fcheck=all -finit-integer=snan -finit-real=snan -finit-logical=true -finit-character=42  -ffpe-trap=invalid,zero,overflow " cmake ../
make -j VERBOSE=1
