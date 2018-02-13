#!/usr/bin/env bash
set -x
rm -rf $HOME/esl
rm -rf $HOME/modules
source /etc/profile.d/lmod.sh



git clone https://gitlab.e-cam2020.eu/esl/esl-bundle.git
pushd esl-bundle 
mkdir -p /home/drFaustroll/esl-bundle/install/include
echo "prefix='/home/drFaustroll/esl/openmpi/0.0.1'" >> rcfiles/debian-gcc-openmpi.rc
sed -i "/pspio/d" rcfiles/debian-gcc-openmpi.rc
echo "module_autogenargs['pspio'] = \"--enable-gsl\"" >> rcfiles/debian-gcc-openmpi.rc
./jhbuild.py -f rcfiles/debian-gcc-openmpi.rc build esl-bundle
popd

mkdir -p $HOME/modules/esl/openmpi
mv modules-esl-openmpi $HOME/modules/esl/openmpi/0.0.1
rm -rf esl-bundle
rm -rf esl/openmpi/0.0.1/_jhbuild

