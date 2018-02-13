#!/usr/bin/env bash
source /etc/profile.d/lmod.sh
module use $HOME/modules
module load gnu/7  fftw3/3.3.6 openmpi/3.0.0 scalapack/2.0.2 esl/openmpi/0.0.1
mkdir build
pushd build
cmake ../ -DWITH_DOC=on
make doc
popd
mkdir -p public/developers
pushd doc/manual
make html
popd
