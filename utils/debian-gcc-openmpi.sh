#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.3.0
mkdir build
cd build
cmake ../
make -j VERBOSE=1
