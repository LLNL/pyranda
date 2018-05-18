################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
################################################################################
import os
import sys
import subprocess
from setuptools import setup
from distutils.command.build import build
from distutils.command.install import install

try:
    import mpi4py
    # needed to test mpi4py+mpif90 in parcop's makefile
    from numpy import f2py
except ImportError as e:
    print(e)
    sys.exit("for now we require you to pre-install mpi4py and numpy")

# try to get the mpif90 compiler from mpi4py, it isn't always there..
# if a user installs like `env MPICC=/path/to/mpicc pip install mpi4py ...`
# it doesn't seem to have anything other than mpicc
mpi4py_compilers = mpi4py.get_config()
if 'mpif90' in mpi4py_compilers:
    mpi4py_mpif90 = mpi4py_compilers['mpif90']
elif 'mpifort' in mpi4py_compilers:
    mpi4py_mpif90 = mpi4py_compilers['mpifort']
# last effort, try to build the possible location
elif 'mpicc' in mpi4py_compilers and os.path.exists(mpi4py_compilers['mpicc'][:-2] + 'f90'):
    mpi4py_mpif90 = mpi4py_compilers['mpicc'][:-2] + 'f90'

distname = "pyranda"
fortran_module = 'parcop'
fortran_package = '{}/{}'.format(distname, fortran_module)
fortran_module_lib = '{}.so'.format(fortran_module)
packages = [distname, fortran_package]
install_requires = ['matplotlib']
description = """
A Python driven, Fortran powered Finite Difference solver for arbitrary
hyperbolic PDE systems. This is the mini-app for the Miranda code.
"""


class PyrandaMakeMixin():
    user_options = [
        ('mpif90=', None, 'mpif90 compiler'),
        ('use-mpi4py-mpif90=', None, 'use mpi4py\'s mpif90 compiler: '.format(mpi4py_mpif90)),
        ('mpiexec=', None, 'mpi exec command used to verify when verifying the mpi compiler'),
        ('numprocs-arg=', None, 'mpi exec num procs arg used when verifying the mpi compiler'),
        ('numprocs=', None, 'number of procs used when verifying mpi'),
        ('no-mpi-compiler-check=', None, 'disable checking that the mpi compiler is compatible with mpi4py')
    ]

    def initialize_options(self):
        self.mpif90 = None
        self.mpiexec = None
        self.numprocs_arg = None
        self.numprocs = None
        self.no_mpi_compiler_check = None
        self.use_mpi4py_mpif90 = None

    def finalize_options(self):
        if self.no_mpi_compiler_check is None:
            self.no_mpi_compiler_check = False
        if self.use_mpi4py_mpif90 is None:
            self.use_mpi4py_mpif90 = False
        if self.use_mpi4py_mpif90 is not None and self.mpif90 is not None:
            raise RuntimeError("cannot define use-mpi4py-mpif90 and mpif90")

    def clean(self):
        print("cleaning up from {} build".format(fortran_module))
        try:
            subprocess.check_call(['make', '-C', fortran_package, 'clean'])
        except subprocess.CalledProcessError:
            print("failed to clean {}".format(fortran_module))
            raise
        print("{} cleaned".format(fortran_module))

    def run(self):
        python = sys.executable
        # build lib*.a
        print("building {}".format(fortran_module))
        try:
            # TODO: we shouldn't build in the source directory
            args = ['make', '-C', fortran_package, 'python={}'.format(python)]

            if self.use_mpi4py_mpif90 is True:
                if mpi4py_mpif90 is not None:
                    args.append('mpif90={}'.format(mpi4py_mpif90))
                else:
                    raise RuntimeError("unable to determine mpi4py's mpif90 compiler")
            if self.mpif90 is not None:
                args.append('mpif90={}'.format(self.mpif90))
            if self.numprocs_arg is not None:
                args.append('np_arg={}'.format(self.numprocs_arg))
            if self.mpiexec is not None:
                args.append('mpirun={}'.format(self.mpiexec))
            if self.numprocs is not None:
                args.append('numprocs={}'.format(self.numprocs))
            if self.no_mpi_compiler_check is True:
                print("skipping mpi4py verification")
                args.append('verified=yes')

            subprocess.check_call(args)
        except subprocess.CalledProcessError:
            print("failed to build {}".format(fortran_module))
            raise
        print("{} built".format(fortran_module))


class BuildPyranda(build, PyrandaMakeMixin):
    user_options = build.user_options + PyrandaMakeMixin.user_options

    def initialize_options(self):
        build.initialize_options(self)
        PyrandaMakeMixin.initialize_options(self)

    def finalize_options(self):
        build.finalize_options(self)
        PyrandaMakeMixin.finalize_options(self)

    def run(self):
        PyrandaMakeMixin.run(self)
        build.run(self)


class InstallPyranda(install, PyrandaMakeMixin):
    def initialize_options(self):
        install.initialize_options(self)
        PyrandaMakeMixin.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)
        PyrandaMakeMixin.finalize_options(self)

    def run(self):
        PyrandaMakeMixin.run(self)
        install.run(self)
        PyrandaMakeMixin.clean(self)


setup_args = dict(
    name=distname,
    description=description,
    packages=packages,
    package_data={fortran_package: [fortran_module_lib]},
    install_requires=install_requires,
    cmdclass={
        'build': BuildPyranda,
        'install': InstallPyranda
    }
)

setup(**setup_args)
