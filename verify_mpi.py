import sys
import argparse

# make sure we have f2py
try:
    from numpy import f2py
except ImportError as e:
    print(e)
    sys.exit("cannot import f2py")


# make sure we can initialize mpi
try:
    from mpi4py import MPI
except ImportError as e:
    print(e)
    sys.exit("cannot import mpi4py")

# init mpi
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

source = b'''
function sum_comm_ranks(comm) result(ranksum)
  use mpi
  implicit none
  integer, intent(in) :: comm
  integer :: ranksum, rank, err
  call MPI_Comm_rank(comm, rank, err)
  call MPI_Allreduce(rank, ranksum, 1, MPI_INT, MPI_SUM, comm, err)
end function sum_comm_ranks
'''

mod_name = 'fort_mod'

parser = argparse.ArgumentParser(description='test mpi4py compatibility')
parser.add_argument('mpif90', help='path to mpif90 compiler')

args = parser.parse_args()

# compile the fort module
status = 0
if rank == 0:
    status = f2py.compile(source,
        extension='.f90',
        extra_args='--f90exec={}'.format(args.mpif90),
        modulename=mod_name)

    if status != 0:
        print("compile failed!")

status = comm.bcast(status, root=0)

if status != 0:
    sys.exit(1)

# try to import the fortran module
try:
    if rank == 0:
        print("trying to import fortran module")
    from fort_mod import sum_comm_ranks
except ImportError as e:
    print(e)
    sys.exit("cannot import fortran module")

# make sure it functions as expected
if rank == 0:
    print("trying to run fortran function")
fsum = sum_comm_ranks(comm.py2f())
expected = sum(range(comm.Get_size()))

if fsum != expected:
    print("rank {}: {} != {}".format(rank, fsum, expected))

assert(fsum == expected)

if rank == 0:
    import glob, shutil, os
    files = glob.glob('{}*'.format(mod_name))
    for f in files:
        print("removing {}".format(f))
        if os.path.isfile(f):
            os.remove(f)
        elif os.path.isdir(f):
            shutil.rmtree(f)
    print("all good!")
