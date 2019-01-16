#!/usr/bin/env bash
set -x
rm -rf $HOME/esl
rm -rf $HOME/modules
source /etc/profile.d/lmod.sh



git clone https://gitlab.e-cam2020.eu/esl/esl-bundle.git
cp -r esl-bundle esl-bundle-mpi
pushd esl-bundle-mpi 
mkdir -p /home/drFaustroll/esl-bundle-mpi/install/include
echo "prefix='/home/drFaustroll/esl/openmpi/0.3.0'" >> rcfiles/debian-gcc-openmpi.rc
sed -i "/'elpa-mpi'/d" rcfiles/debian-gcc-openmpi.rc
echo "module_autogenargs['elpa-mpi']=autogenargs+\" CFLAGS='-march=native' FCFLAGS='-ffree-line-length-none' --disable-avx2 --disable-avx --disable-sse-assembly \" " >> rcfiles/debian-gcc-openmpi.rc 
./jhbuild.py -f rcfiles/debian-gcc-openmpi.rc build esl-bundle-mpi
popd

pushd esl-bundle
echo "prefix='/home/drFaustroll/esl/serial/0.3.0'" >> rcfiles/debian-gcc-serial.rc
sed -i "/'elpa'/d" rcfiles/debian-gcc-serial.rc
echo "module_autogenargs['elpa']=autogenargs+\" CFLAGS='-march=native' FCFLAGS='-ffree-line-length-none' --disable-avx2 --disable-avx --disable-sse-assembly \" " >> rcfiles/debian-gcc-serial.rc 
./jhbuild.py -f rcfiles/debian-gcc-serial.rc build esl-bundle
popd

mkdir -p $HOME/modules/esl/openmpi
mv modules-esl-openmpi $HOME/modules/esl/openmpi/0.3.0
rm -rf esl-bundle-mpi
rm -rf esl/openmpi/0.3.0/_jhbuild

mkdir -p $HOME/modules/esl/serial
mv modules-esl-serial $HOME/modules/esl/serial/0.3.0
rm -rf esl-bundle
rm -rf esl/serial/0.3.0/_jhbuild

