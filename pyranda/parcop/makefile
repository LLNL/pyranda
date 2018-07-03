#===================================================================================================
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
#===================================================================================================
#
# It's important to check that mpi4py is using the same mpi as what's used here,
# there is a test-mpi4py target that will, given a python and mpif90, try to build and load a
# fortran function from python to verify that the combination is okay.
# You can run `python -c "import mpi4py;print mpi4py.get_config()"` to see which compilers
# mpi4py is using, which will let you get at the underlying mpi.
#
# parallel make does not work because of modfile dependencies in the f90 files...
#
# TODO: it seems like f2py is too fragile for systems and/or compilers that aren't very normal
# so itstead f2py is just used to generate the wrapper source and we define how to compile
# that. A config/host file is probably the way to go.
#
# it would need:
# - mpif90
# - cc
# - fflags
#  * position independent code
#  * free form fortran
# - c flags
#  * position independent code
# - link flags
#  * shared library option (usually -shared)
#
# XXX: anything else??
#

ifeq ($(SYS_TYPE),)
  systype = $(shell uname)
else
	systype = $(SYS_TYPE)
endif

chaos_5_x86_64_ib.gnu.mpif90   = /collab/usr/gapps/python/build/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.4.7/mpich-3.0.4-evwtvwqv5i2ujysfnivcdtrnppmpcsrg/bin/mpif90
chaos_5_x86_64_ib.intel.mpif90 = /usr/local/tools/mvapich2-intel-2.2/bin/mpif90
chaos_5_x86_64_ib.python       = /usr/local/bin/python
chaos_5_x86_64_ib.compiler     = intel
chaos_5_x86_64_ib.mpirun       = srun
chaos_5_x86_64_ib.np_arg       = -n

toss_3_x86_64_ib.gnu.mpif90   = /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-6.1.0/bin/mpif90
toss_3_x86_64_ib.intel.mpif90 = /usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.3/bin/mpif90
toss_3_x86_64_ib.clang.mpif90 = /usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpif90
toss_3_x86_64_ib.pgi.mpif90   = /usr/tce/packages/mvapich2/mvapich2-2.2-pgi-18.1/bin/mpif90
toss_3_x86_64_ib.python       = /usr/tce/bin/python
toss_3_x86_64_ib.compiler     = intel
toss_3_x86_64_ib.mpirun       = srun
toss_3_x86_64_ib.np_arg       = -n

blueos_3_ppc64le_ib.xl.mpif90  = /usr/tce/packages/spectrum-mpi/spectrum-mpi-2018.04.27-xl-beta-2018.06.01/bin/mpif90
blueos_3_ppc64le_ib.xl.cc      = /usr/tce/packages/xl/xl-beta-2018.06.01/bin/xlc
blueos_3_ppc64le_ib.gnu.mpif90 = /usr/tcetmp/bin/mpigfortran
blueos_3_ppc64le_ib.compiler   = gnu
blueos_3_ppc64le_ib.python     = /usr/tcetmp/bin/python
blueos_3_ppc64le_ib.mpirun     = mpirun
blueos_3_ppc64le_ib.np_arg     = -n

Darwin.gnu.mpif90 = $(shell which mpif90)
Darwin.python     = $(shell which python)
Darwin.compiler   = gnu
Darwin.mpirun     = mpirun
Darwin.np_arg     = -n

prom.pgi.mpif90 = $(shell which mpif90)
prom.python     = $(shell which python)
prom.compiler   = pgi
prom.mpirun     = mpirun
prom.np_arg     = -n
prom.pgi.cc     = pgcc
prom.gnu.cc     = g++

# if the systype.compiler is empty setup will use the following:
default.gnu.mpif90 = $(shell which mpif90)
default.python     = $(shell which python)
default.compiler   = gnu
default.mpirun     = $(shell which mpirun)
default.np_arg     = -n

gnu.fflags   = -fPIC -ffree-form -ffree-line-length-0 -fbackslash
intel.fflags = -fPIC
# TODO: what are clang's fflags here
clang.fflags = -fPIC
xl.fflags    = -fPIC -qsuppress=cmpmsg -qsuppress=1500-036 -qxlf2003=polymorphic -qnosave -qfree=f90 -qsuffix=cpp=f:f=f -qmaxmem=-1 -qextname
xl.lopts     = -qmkshrobj
pgi.fflags   = -fPIC
ifdef gpu
pgi.fflags += CUDAROOT=/gpfs/sw/software/CUDA/9.1.85/ -Mcuda -ta=tesla:cc70,managed -Minfo=accel
pgi.ldflags += CUDAROOT=/gpfs/sw/software/CUDA/9.1.85/ -Mcuda -ta=tesla:cc70,managed -Minfo=accel
pgi.extlibs += -lcublas
endif

default.lopts = -shared

ifeq ($($(systype).compiler),)
	systype=default
endif

compiler ?= $($(systype).compiler)
mpirun   ?= $($(systype).mpirun)
np_arg   ?= $($(systype).np_arg)
nprocs   ?= 2
python   ?= $($(systype).python)
cc       ?= $(if $($(systype).$(compiler).cc),$($(systype).$(compiler).cc),cc)
mpif90   ?= $($(systype).$(compiler).mpif90)
# make sure we use this python's f2py unless the user wants something else
f2py     ?= $(python) -c 'import numpy.f2py; numpy.f2py.main()'
fflags   ?= $($(compiler).fflags)
lopts    ?= $(if $($(compiler).lopts),$($(compiler).lopts),$(default.lopts))
ldflags  ?= $($(compiler).ldflags)
extlibs  ?= $($(compiler).extlibs)

ifeq ($(profile),1)
fflags += -pg
ldflags += -pg
endif
#
# The python/mpif90 combination is saved to a file (if the mpi4py integration test passes)
# so that we don't have to test everytime a user calls make.
# the parcop.so target only uses the test-mpi4py target if the py/mpi combo hasn't 
# been verified, otherwise it always rebuilds because test-mpi4py is never "up-to-date".
#
verifile  = verified_compiler
verinfo   = "$(python),$(mpif90)"
verified ?= $(if $(shell test -f $(verifile) && \
					 test "`cat $(verifile)`" = "$(verinfo)" && echo "yes")\
					 ,yes,no)
ifneq ($(verified),yes)
	verified_target=test-mpi4py
endif

no_mpirun_verify_cmd = $(python) verify_mpi.py $(mpif90)
mpirun_verify_cmd    = $(mpirun) $(np_arg) $(nprocs) $(no_mpirun_verify_cmd)
verify_cmd 				   = $(mpirun_verify_cmd)

#
# Python and numpy directories for use when building the f2py wrapped module
# TODO: sysconfig.get_config_vars() has a lot of stuff that might make this better
#
python.Iflags = -I$(shell $(python) -c "from sysconfig import get_config_var; print(get_config_var('INCLUDEPY'))")
python.Lflags = -L$(shell $(python) -c "from sysconfig import get_config_var; print(get_config_var('LIBDIR'))")
python.major  = $(shell $(python) -c "import sys; print(sys.version_info[0])")
python.minor  = $(shell $(python) -c "import sys; print(sys.version_info[1])")
# TODO: a better way (use one of the vars in sysconfig.get_config_vars instead)
# it's libpython2.7 or libpython3, libpython3.6m
ifeq ($(python.major),2)
	python.lflags = -lpython$(python.major).$(python.minor)
else
	python.lflags = -lpython$(python.major)
endif
#numpy.dir    = $(shell $(python) -c "from sysconfig import get_paths; print(get_paths()['purelib'])")/numpy
numpy.dir    = $(shell $(python) -c "import numpy; print(numpy.__path__[0])")
# TODO: python -c "import numpy; print(numpy.__path__)" might be better
numpy.Iflags = -I$(numpy.dir)/core/include
f2py.dir     = $(numpy.dir)/f2py/src
f2py.Iflags  =-I$(f2py.dir)

#
# source order is important because of mod file dependencies
#
source = blockmath.o stencils.o patch.o pentadiagonal.o \
	     comm.o compact.o mesh.o objects.o compact_operators.o \
			 operators.o

vpath %.c $(f2py.dir)

.PHONY: clean info test-mpi4py

%.o: %.f90
	$(mpif90) $(fflags) -c $<

%.o: %.c
	echo $(cc) -fPIC $(python.Iflags) $(numpy.Iflags) $(f2py.Iflags) -c $^
	$(cc) -fPIC $(python.Iflags) $(numpy.Iflags) $(f2py.Iflags) -c $^

all: parcop.so test_pent test_laplace

# See nasa's page for reference:
#	 https://modelingguru.nasa.gov/docs/DOC-2412
#parcop.so: libparcop.a parcop.f90
#	@$(f2py) \
#		-c --f90exec=$(mpif90) \
#		-m parcop \
#		parcop.f90

libparcop.a: $(source)
	@ar -ur $@ $(source)

# TODO: is the '2' in of the f2pywrappers2.f90 the major python version?
#                                        ^
parcopmodule.c parcop-f2pywrappers2.f90: libparcop.a parcop.o
	@$(f2py) --lower --wrap-functions -m parcop parcop.f90 --verbose

parcop.so: parcopmodule.o fortranobject.o parcop-f2pywrappers2.o
	$(mpif90) $(lopts) $(python.Lflags) \
		parcopmodule.o fortranobject.o parcop.o parcop-f2pywrappers2.o libparcop.a $(python.lflags) \
		-o parcop.so


#
# Checks to see that the python,mpif90 combo has been verified, if not
# try to verify. If the verification is successful write the combo to the
# "verifile" otherwise fail
#
test-mpi4py:
	@if [ "$(verified)" = "yes" ]; then \
		echo "compiler already verified"; \
	else \
		echo "running: $(verify_cmd)"; \
		$(verify_cmd) && { echo "$(verinfo)" > $(verifile); echo "$(mpif90) verified"; } || \
			{ echo "$(mpif90) failed mpi4py verification"; false; }; \
	fi

test_pent: libparcop.a parcop.o test_pent.o
	$(mpif90) $(ldflags) test_pent.o parcop.o libparcop.a -o $@ $(extlibs)

test_laplace: libparcop.a parcop.o test_laplace.o
	$(mpif90) $(ldflags) test_laplace.o parcop.o libparcop.a -o $@ $(extlibs)

info:
	@echo "systype:    $(systype)"
	@echo "python:     $(python) - verified: $(verified)"
	@echo "mpif90:     $(mpif90) - verified: $(verified)"
	@echo "f2py:       $(f2py)"
	@echo "fflags:     $(fflags)"
	@echo "mpirun:     $(mpirun)"
	@echo "nprocs arg: $(np_arg)"
	@echo "nprocs:     $(nprocs)"
	@echo "compiler:   $(compiler)"

clean:
	rm -rf *.o *.mod test_pent *.so *.pyf *.pyc libparcop.a *.dSYM $(verifile) parcopmodule.c parcop-f2pywrappers*.f90
