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
# It's important to check that mpi4py is using the same mpi as what you build with here
# You can run `python -c "import mpi4py;print mpi4py.get_config()"`, which will confirm that
# on toss3 machines mpi4py is built with mvapich2-2.2.
#
# parallel make does not work because of modfile dependencies in the f90 files...
#
# no -fcompiler option is given to f2py although, at least with intel, you can
# pass this. I'm not sure if it's important or not, that's something to look at later.
#


systype = $(SYS_TYPE)


ifeq ($(systype), )
  systype = $(shell uname)
endif


ifeq ($(systype),chaos_5_x86_64_ib)


  cpath = /collab/usr/gapps/python/build/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.4.7/mpich-3.0.4-evwtvwqv5i2ujysfnivcdtrnppmpcsrg/bin
  gnu.f77 = $(cpath)/mpif77
  gnu.f90 = $(cpath)/mpif90
  gnu.fflags = -ffree-form -ffree-line-length-0 -fbackslash

  ipath = /usr/local/tools/mvapich2-intel-2.2
  intel.f77 = $(ipath)/bin/mpif77
  intel.f90 = $(ipath)/bin/mpif90

  compiler = intel

  chaos_5_x86_64_ib.f2py = f2py

else ifeq ($(systype),toss_3_x86_64_ib)

  gnu.f77 = /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-6.1.0/bin/mpif77
  gnu.f90 = /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-6.1.0/bin/mpif90
  gnu.fflags = -ffree-form -ffree-line-length-0 -fbackslash

  intel.f77 = /usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.3/bin/mpif77
  intel.f90 = /usr/tce/packages/mvapich2/mvapich2-2.2-intel-16.0.3/bin/mpif90

  clang.f77 = /usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpif77
  clang.f90 = /usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpif90
  clang.fflags = -ffree-form -ffree-line-length-0 -fbackslash

  compiler = intel

  toss3_x86_64_ib.f2py = f2py


else ifeq ($(systype),Darwin)

  gnu.f77 = /usr/local/bin/mpif77
  gnu.f90 = /usr/local/bin/mpif90

  comp := $(shell command -v $(gnu.f90))
  ifndef comp
    gnu.f77 = /opt/local/bin/mpif77
    gnu.f90 = /opt/local/bin/mpif90
  endif


  gnu.fflags = -ffree-form -ffree-line-length-0 -fbackslash	

  compiler = gnu

  Darwin.f2py = f2py-2.7


endif 

f90 ?= $($(compiler).f90)
f2py ?= $($(systype).f2py)

# Position independent code 
fflags = -fPIC $($(compiler).fflags)
ifeq ($(DEBUG),1)
  fflags = -fPIC $($(compiler).fflags) -g -O0
endif

exec = PENT

source = blockmath.o stencils.o patch.o pentadiagonal.o \
	     comm.o compact.o mesh.o objects.o compact_operators.o operators.o parcop.o

.PHONY: clean

%.o: %.f90
	$(f90) $(fflags) -c $<

# See nasa's page for reference:
#	 https://modelingguru.nasa.gov/docs/DOC-2412
all: libparcop.a fort-test
	$(f2py) \
		--include-paths . \
		-c --f90exec=$(f90) \
		libparcop.a \
		-m parcop \
		parcop.f90

libparcop.a: $(source)
	ar -ur $@ $(source) 

fort-test : test_pent.o libparcop.a
	$(f90) $^ $(fflags) -o $@

check-mpi4py:
	@echo $(systype)
	python -c "import mpi4py;print mpi4py.get_config()"


clean:
	rm -f *.o *.mod $(exec) *.so *.pyf *.pyc libparcop.a fort-test
