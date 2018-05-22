# Installing Pyranda

## Prerequisites
Before installing pyranda you'll need to install numpy and mpi4py. Numpy's f2py utility
is used to build the python-fortran interface and we try to validate that your choice of
compiler will work with the installed mpi4py's version of mpi.

### installing `mpi4py`
`mpi4py` should be installed with the `--no-cache-dir` option to avoid using an
existing build with a cached compiler.

```
# if your mpi*s are in your path
pip install mpi4py --no-cache-dir
# otherwise you can specify an environment variable
env MPICC=/path/to/your/mpi pip install mpi4py --no-cache-dir
```

### installing `numpy`
`numpy` shouldn't have any special build steps, just install as normal:

```
pip install numpy
```

## Installing with pip

### [optional] Using a virtualenv

```
...> virtualenv [-p python] my_venv
...> source my_venv/bin/activate
(my_venv) ...>
```

You can also verify that you're in your venv by checking your `$PATH`:

```
...> echo $PATH
/path/to/your/env/root/my_venv/bin:...
```

### install pyranda

```
pip install . [--user]
```

## Installing without pip

```
python setup.py install
```


## Legacy Instructions


### python-
Though other versions of python may very well work, we recommend and support
python 2.7 for pyranda.  Python will need to see the python of pyranda.

** Add pyranda/src/python to PYTHONPATH **

### numpy-
As long as numpy is working with your version of python above, there will be no
compability issues.  This can be installed in a number of ways. http://www.numpy.org

### mpi4py
This python package provides MPI bindings to python and may or may not exists on your system
and python path.  At a minimum, the version of MPI used to build mpi4py must match the version
of MPI of your fortran compiler.  

To check your mpi4py config, run:
python -c "import mpi4py;print mpi4py.get_config()"

It is likely that you will need to build your own version of mpi4py to suite your specified compiler
type.  https://bitbucket.org/mpi4py/mpi4py

### Example install (this should work on LLNL-LC)
`wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz`
`tar xvzf mpi4py-3.0.0.tar.gz`
`cd mpi4py`

#### Add this to mpi.cfg file
[llnl-intel]
mpi_dir              = /usr/local/tools/mvapich2-intel-2.2
mpicc                = %(mpi_dir)s/bin/mpicc
mpicxx               = %(mpi_dir)s/bin/mpicxx
library_dirs         = %(mpi_dir)s/lib


`python setup.py build --mpi=llnl-intel`
`python setup.py install --prefix=install_location`

** Add install_location/*/site_packages to PYTHONPATH **


### fortran-
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
