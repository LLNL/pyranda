################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
################################################################################
from mpi4py import MPI
import sys
import argparse
import glob
import shutil
import os
from numpy import f2py

# init mpi
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# the fortran function we try to run
source = b'''
function sum_comm_ranks(comm) result(ranksum)
  use mpi
  implicit none
  integer, intent(in) :: comm
  integer :: ranksum, mpirank, err
  call MPI_Comm_rank(comm, mpirank, err)
  call MPI_Allreduce(mpirank, ranksum, 1, MPI_INT, MPI_SUM, comm, err)
end function sum_comm_ranks
'''

# parse command line
status = 0
if rank == 0:
    parser = argparse.ArgumentParser(description='test mpi4py compatibility')
    parser.add_argument('mpif90', help='path to mpif90 compiler')

    try:
        args = parser.parse_args()
    except SystemExit:
        status = 1

if comm.bcast(status, root=0) != 0:
    sys.exit(1)


# compile the fort module
status = 0
# mod_name is used below in an import
mod_name = 'fort_mod'
if rank == 0:
    status = f2py.compile(source,
        extension='.f90',
        extra_args='--f90exec={}'.format(args.mpif90),
        verbose=False,
        modulename=mod_name)

    if status != 0:
        print("compile failed!")

if comm.bcast(status, root=0) != 0:
    sys.exit(1)


def print0(arg):
    if rank == 0:
        print(arg)


def cleanup():
    if rank == 0:
        files = glob.glob('{}*'.format(mod_name))
        for f in files:
            print("removing {}".format(f))
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)


# try to import the fortran module
try:
    print0("trying to import fortran module...")
    from fort_mod import sum_comm_ranks
except ImportError as e:
    print0(e)
    sys.exit("cannot import fortran module")
finally:
    cleanup()

print0("successfully imported")


# make sure it functions as expected
print0("trying to run fortran function...")
# tries to make the output ordering nice
comm.barrier()
fsum = sum_comm_ranks(comm.py2f())
expected = sum(range(comm.Get_size()))
print("rank {}: expected: {}, got: {}".format(rank, fsum, expected))
if fsum != expected:
    print("expected != got")
else:
    print0("successfully ran")

assert(fsum == expected)

if rank == 0:
    print("all good!")
