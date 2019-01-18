#!/usr/bin/env bash
set -x
rm -rf $HOME/esl
rm -rf $HOME/modules
source /etc/profile.d/lmod.sh

module load gnu/8 openmpi/3.1.1 scalapack/2.0.2 pfftw3/3.3.6
export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH
git clone https://gitlab.com/ElectronicStructureLibrary/esl-bundle.git
cp -r esl-bundle esl-bundle-mpi
pushd esl-bundle-mpi
mkdir -p /home/drFaustroll/esl-bundle-mpi/install/include
echo "prefix='/home/drFaustroll/esl/openmpi/0.3.1'" >> rcfiles/opensuse-gcc-openmpi.rc
sed -i "/'elpa-mpi'/d" rcfiles/opensuse-gcc-openmpi.rc
echo "module_autogenargs['elpa-mpi']=autogenargs+\" CFLAGS='-march=native' FCFLAGS='-ffree-line-length-none' --disable-avx2 --disable-avx --disable-sse-assembly \" " >> rcfiles/opensuse-gcc-openmpi.rc 
./jhbuild.py -f rcfiles/opensuse-gcc-openmpi.rc build esl-bundle-mpi
popd
module purge
unset LIBRARY_PATH
module load gnu/8 fftw3/3.3.6
pushd esl-bundle
mkdir -p /home/drFaustroll/esl-bundle/install/include
echo "prefix='/home/drFaustroll/esl/serial/0.3.1'" >> rcfiles/opensuse-gcc-serial.rc

sed -i "/'elpa'/d" rcfiles/opensuse-gcc-serial.rc
echo "module_autogenargs['elpa']=autogenargs+\" CFLAGS='-march=native' FCFLAGS='-ffree-line-length-none' --disable-avx2 --disable-avx --disable-sse-assembly \" " >> rcfiles/opensuse-gcc-serial.rc 

./jhbuild.py -f rcfiles/opensuse-gcc-serial.rc build esl-bundle
popd

mkdir -p $HOME/modules/esl/openmpi
mv modules-esl-openmpi $HOME/modules/esl/openmpi/0.3.1
rm -rf esl-bundle-mpi
rm -rf esl/openmpi/0.3.1/_jhbuild

mkdir -p $HOME/modules/esl/serial
mv modules-esl-serial $HOME/modules/esl/serial/0.3.1
rm -rf esl-bundle
rm -rf esl/serial/0.3.1/_jhbuild

