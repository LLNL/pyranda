#!/bin/sh

here=$(pwd)
install_dir=${here}/${1}
test -e ${install_dir}/lib/libmpi.so && { echo "$1 already installed"; exit 0; }

case "$1" in
  mpich)
    name=mpich-3.2.1
    wget http://www.mpich.org/static/downloads/3.2.1/${name}.tar.gz
    tar -xzf ${name}.tar.gz
    cd $name
    ./configure --prefix=$install_dir
    make -j4
    make install
    rm -fr ${name}*
  ;;

  openmpi)
    name=openmpi-3.1.0
    wget https://download.open-mpi.org/release/open-mpi/v3.1/${name}.tar.gz
    tar -xzf ${name}.tar.gz
    cd $name
    ./configure --prefix=$install_dir
    make -j4
    make install
    rm -rf ${name}*
  ;;

  *)
    echo "unknown mpi: '$1'"
    exit 1
  ;;
esac
