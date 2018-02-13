#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.0.1
mkdir build
cd build
FFLAGS="-g -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=snan -finit-real=snan -finit-logical=true -finit-character=42 -finit-derived -ffpe-trap=invalid,zero,overflow -fdump-core
" cmake ../
make -j VERBOSE=1
