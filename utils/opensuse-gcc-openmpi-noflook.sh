#!/usr/bin/env bash
source /etc/profile.d/lmod.sh
module use $HOME/modules
module load gnu/7  fftw3/3.3.6 openmpi/3.0.0 scalapack/2.0.2 esl/openmpi/0.1.0
export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH
mkdir build
pushd build
cmake ../ -DWITH_FLOOK=Off
make -j VERBOSE=1
