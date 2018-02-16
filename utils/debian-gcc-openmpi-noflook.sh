#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.1.0
mkdir build
cd build
cmake ../ -DWITH_FLOOK=Off
make -j VERBOSE=1
