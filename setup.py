import os
import sys
import subprocess
from setuptools import setup
from distutils.command.build import build
from distutils.command.install import install

try:
    import mpi4py
    from numpy import f2py
except ImportError as e:
    print(e)
    sys.exit("for now we require you to pre-install mpi4py and numpy")

distname = "pyranda"
fortran_module = 'parcop'
fortran_package = '{}/{}'.format(distname, fortran_module)
fortran_module_target = 'libparcop.a'
fortran_module_interface = 'parcop.f90'
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
        ('mpiexec=', None, 'mpi launcher command'),
        ('numprocs-arg=', None, 'mpi launcher nprocs arg')
    ]

    def initialize_options(self):
        self.mpif90 = None
        self.mpiexec = None
        self.numprocs_arg = None

    def finalize_options(self):
        if self.mpif90 is None:
            self.mpif90 = 'mpif90'
        if self.mpiexec is None:
            self.mpiexec = 'mpirun'
        if self.numprocs_arg is None:
            self.numprocs_arg = '-n'

    def run(self):
        python = sys.executable

        print("verifying that mpi4py is compatiable with your mpifortran compiler")
        try:
            subprocess.check_call([self.mpiexec, self.numprocs_arg, '4', python, 'verify_mpi.py', self.mpif90])
        except subprocess.CalledProcessError:
            print("failed to verify mpi")
            raise
        print("{} is compatible".format(self.mpif90))

        print("building {}".format(fortran_module_target))
        try:
            subprocess.check_call(['make', '-C', fortran_package, 'f90={}'.format(self.mpif90), 'f2py=f2py'])#fortran_module_target])
        except subprocess.CalledProcessError:
            print("failed to build {}".format(fortran_module_target))
            raise
    
        '''
        print("building parcop's python interface")
        with open(os.path.join(fortran_package, fortran_module_interface), 'r') as _:
            source = _.read()
        status = f2py.compile(source,
            extension='.f90',
            extra_args='--f90exec={} {} --include-paths ./pyranda/parcop'.format(self.mpif90, fortran_module_target),
            modulename=fortran_module)
        if status != 0:
            raise RuntimeError('failed to compile the fortran python interface')
        '''


class BuildPyranda(build, PyrandaMakeMixin):
    user_options = build.user_options + PyrandaMakeMixin.user_options

    def initialize_options(self):
        build.initialize_options(self)
        PyrandaMakeMixin.initialize_options(self)

    def finalize_options(self):
        build.finalize_options(self)
        PyrandaMakeMixin.finalize_options(self)

    def run(self):
        build.run(self)
        PyrandaMakeMixin.run(self)


class InstallPyranda(install, PyrandaMakeMixin):
    user_options = install.user_options + PyrandaMakeMixin.user_options

    def initialize_options(self):
        install.initialize_options(self)
        PyrandaMakeMixin.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)
        PyrandaMakeMixin.finalize_options(self)

    def run(self):
        install.run(self)
        PyrandaMakeMixin.run(self)


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
