FROM debian:latest

RUN groupadd -g 999 drFaustroll && useradd -r -u 999 -g drFaustroll drFaustroll && mkdir /home/drFaustroll && chown -R drFaustroll:drFaustroll \ 
/home/drFaustroll && apt update &&  apt install -y libscalapack-mpi-dev cmake make autoconf automake libtool \
python sudo git wget curl libblas-dev liblapack-dev libhdf5-openmpi-dev pkg-config liblua5.2-dev \
lmod procps check doxygen libhdf5-serial-dev libfftw3-dev libgsl-dev g++ vim nano
WORKDIR /home/drFaustroll
USER drFaustroll
COPY assets/*  /home/drFaustroll/
RUN ./newsoftware.sh && rm -f ./newsoftware.sh && rm -f ./modules-* 
