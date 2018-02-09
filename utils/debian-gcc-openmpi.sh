#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.0.1
mkdir build
cd build
cmake ../
make -j
