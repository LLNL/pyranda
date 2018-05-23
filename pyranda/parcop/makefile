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
# TODO: additional f2py flags we might want to look at:
#  --f77flags=          Specify F77 compiler flags
#  --f90flags=          Specify F90 compiler flags
#  --opt=               Specify optimization flags
#  --arch=              Specify architecture specific optimization flags
#  --noopt              Compile without optimization
#  --noarch             Compile without arch-dependent optimization
#  --debug              Compile with debugging information
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
toss_3_x86_64_ib.python       = /usr/tce/bin/python3
toss_3_x86_64_ib.compiler     = intel
toss_3_x86_64_ib.mpirun       = srun
toss_3_x86_64_ib.np_arg       = -n

blueos_3_ppc64le_ib.ibm.mpif90 = /usr/tcetmp/packages/spectrum-mpi/spectrum-mpi-2017.04.03/bin/mpif90
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
xlf.fflags = -fPIC -qsuppress=cmpmsg -qsuppress=1500-036 -qxlf2003=polymorphic -qnosave -qfree=f90 -qsuffix=cpp=f:f=f -qmaxmem=-1

ifeq ($($(systype).compiler),)
	systype=default
endif

compiler ?= $($(systype).compiler)
mpirun   ?= $($(systype).mpirun)
np_arg   ?= $($(systype).np_arg)
nprocs   ?= 2
python   ?= $($(systype).python)
mpif90   ?= $($(systype).$($(systype).compiler).mpif90)
# make sure we use this python's f2py unless the user wants something else
f2py     ?= $(python) -c 'import numpy.f2py; numpy.f2py.main()'
fflags   ?= $(if $($(systype).$($(systype).$(compiler)).fflags)\
              ,$($(systype).$($(systype).$(compiler)).fflags)\
              ,$($(compiler).fflags))

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

exec = fort-test
source = blockmath.o stencils.o patch.o pentadiagonal.o \
	     comm.o compact.o mesh.o objects.o compact_operators.o operators.o parcop.o

.PHONY: clean info test-mpi4py

%.o: %.f90
	$(mpif90) $(fflags) -c $<

all: parcop.so fort-test

# See nasa's page for reference:
#	 https://modelingguru.nasa.gov/docs/DOC-2412
parcop.so: $(verified_target) libparcop.a parcop.f90
	@$(f2py) \
		--include-paths . \
		-c --f90exec=$(mpif90) \
		libparcop.a \
		-m parcop \
		parcop.f90

libparcop.a: $(source)
	@ar -ur $@ $(source)

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

$(exec): test_pent.o libparcop.a
	$(mpif90) $^ $(fflags) -o $@

info:
	@echo "systype:    $(systype)"
	@echo "python:     $(python) - verified: $(verified)"
	@echo "mpif90:     $(mpif90) - verified: $(verified)"
	@echo "f2py:       $(f2py)"
	@echo "fflags:     $(fflags)"
	@echo "mpirun:     $(mpirun)"
	@echo "nprocs arg: $(np_arg)"
	@echo "nprocs:     $(nprocs)"

clean:
	rm -rf *.o *.mod $(exec) *.so *.pyf *.pyc libparcop.a fort-test *.dSYM $(verifile)
