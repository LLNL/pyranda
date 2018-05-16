import os
import sys
from setuptools import setup

try:
    import mpi4py
    import numpy
except ImportError as e:
    print(e)
    sys.exit("for now we require you to pre-install mpi4py and numpy")

parcop_lib = 'pyranda/parcop/parcop.so'
if not os.path.exists(parcop_lib):
    sys.exit("for now we require you to pre-build parcop")

distname = "pyranda"
packages = [distname, '{}/parcop'.format(distname)]
install_requires = ['matplotlib']

setup_args = dict(
    name=distname,
    description="Pyranda",
    packages=packages,
    package_data={'pyranda/parcop': ['parcop.so']},
    install_requires=install_requires
)

setup(**setup_args)
