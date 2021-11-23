#!/bin/bash

install_dir=$2
test -e ${install_dir}/lib/libmpi.so && { echo "$1 already installed"; exit 0; }

case "$1" in
  mpich-*)
    # https://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
    # https://www.mpich.org/static/downloads/3.4/mpich-3.4.tar.gz
    ver="$(echo $1 | cut -d- -f2)"
    wget http://www.mpich.org/static/downloads/${ver}/${1}.tar.gz
    tar -xzf ${1}.tar.gz
    cd $1
    ./configure --prefix=$install_dir --with-device=ch3
    make -j4
    make install
  ;;

  openmpi-*)
    # https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
    ver="$(echo $1 | cut -d- -f2 | cut -d. -f1,2)"
    wget https://download.open-mpi.org/release/open-mpi/v${ver}/${1}.tar.gz
    tar -xzf ${1}.tar.gz
    cd $1
    ./configure --prefix=$install_dir
    make -j4
    make install
  ;;

  *)
    echo "unknown mpi: '$1'"
    exit 1
  ;;
esac