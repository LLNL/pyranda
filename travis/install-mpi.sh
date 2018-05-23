#!/bin/sh

case "$1" in
  mpich)
    name=mpich-3.2.1
    wget http://www.mpich.org/static/downloads/3.2.1/${name}.tar.gz
    tar -xzf ${name}.tar.gz
    cd $name
    ./configure --prefix=$(pwd)/../../mpich
    make -j4
    make install
    rm -fr ${name}*
  ;;

  openmpi)
    name=openmpi-3.1.0
    wget https://download.open-mpi.org/release/open-mpi/v3.1/${name}.tar.gz
    tar -xzf ${name}.tar.gz
    cd $name
    ./configure --prefix=$(pwd)/../../openmpi
    make -j4
    make install
    rm -rf ${name}*
  ;;

  *)
    echo "unknown mpi: '$1'"
    exit 1
  ;;
esac
