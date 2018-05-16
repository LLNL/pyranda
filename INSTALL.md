# Install notes for pyranda

The following outlines what/how to install pyranda.  At a minimum, you will need to 
build one fortran library, and add two paths to your PYTHONPATH.

## Installing with pip
For now you're required to install some python packages and build libparcop
before installing. You can install mpi4py and numpy with pip and it's important
to remember that you must use the same mpi for mpi4py and libparcop. If you want
to use a non-default mpi for mpi4py see the instructions below.

### Building libparcop
There are some default hostname-based mpif90 and f2py paths in *pyranda/parcop/makefile* but
it's probably best to provide your own for now.

```shell
make -C pyranda/parcop [f90=] [f2py=]
```

### Installing with pip
*pip will fail if the prereqs don't exist*

```shell
pip install . [--user]
```

### Running the tests

```
cd tests
python run_tests.py

...

heat1D-analytic-256 -- 2.44249065418e-15



=====Testing Summary======
Passed: 20
Failed: 0



===== New baselines =====
```

## python-
Though other versions of python may very well work, we recommend and support
python 2.7 for pyranda.  Python will need to see the python of pyranda.

** Add pyranda/src/python to PYTHONPATH **

## numpy-
As long as numpy is working with your version of python above, there will be no
compability issues.  This can be installed in a number of ways. http://www.numpy.org

## mpi4py
This python package provides MPI bindings to python and may or may not exists on your system
and python path.  At a minimum, the version of MPI used to build mpi4py must match the version
of MPI of your fortran compiler.  

To check your mpi4py config, run:
python -c "import mpi4py;print mpi4py.get_config()"

It is likely that you will need to build your own version of mpi4py to suite your specified compiler
type.  https://bitbucket.org/mpi4py/mpi4py

## Example install (this should work on LLNL-LC)
`wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz`
`tar xvzf mpi4py-3.0.0.tar.gz`
`cd mpi4py`

### Add this to mpi.cfg file
[llnl-intel]
mpi_dir              = /usr/local/tools/mvapich2-intel-2.2
mpicc                = %(mpi_dir)s/bin/mpicc
mpicxx               = %(mpi_dir)s/bin/mpicxx
library_dirs         = %(mpi_dir)s/lib


`python setup.py build --mpi=llnl-intel`
`python setup.py install --prefix=install_location`

** Add install_location/*/site_packages to PYTHONPATH **


## fortran-
A fortran compiler with 2003 and above standards enforced and MPI libraries is required.
The MPI version used here should match that of mpi4py.

pyranda/src/fortran contains the fortran files and a makefile to create the python extension.  
The fortran compiler will need to be set here.  

f2py is also invoked at build time.  This utility installs with numpy and ought to be in your path.

To build, type:

`make`

Optionally, build debug via:
`make DEBUG=1`

make will produce:
- libparcop.a (fortran library)
- parcop.so (python library which wraps libparcop.a)
- fort-test (fortran exec for simple tests, not used by pyranda)

** Add pyranda/src/fortran to PYTHONPATH **


## Setup guide
To facilitate setup, a setup script will be created in due time which will:

- check for existance of python
- check for existance of numpy, f2py
- check for existance of fortran compilers
   - prompt for selection
- check for mpi4py version
  - if not matched with desired fortran, then will build this automatically
- build pyranda
- set PYTHONPATH 

