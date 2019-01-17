#!/usr/bin/env bash
source /etc/profile.d/lmod.sh
module use $HOME/modules
module load gnu/8  fftw3/3.3.6 openmpi/3.1.1 scalapack/2.0.2 esl/openmpi/0.3.1
export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH
mkdir build
pushd build
cmake ../
make -j VERBOSE=1
