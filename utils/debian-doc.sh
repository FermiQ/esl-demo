#!/usr/bin/env bash

. /etc/profile.d/lmod.sh
module use $HOME/modules
module load esl/openmpi/0.3.1
mkdir build
cd build
cmake ../ -DWITH_DOC=On
make -j
make doc

