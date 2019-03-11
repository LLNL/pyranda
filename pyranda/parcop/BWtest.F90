!pgf90 -o BWtest -ta=tesla,cuda9.2,lineinfo,keepgpu,safecache BWtest.f90 -Minfo=accel -Mpreprocess
!mpixlf -qsmp=omp -qoffload -qfree=f90 -o BWtest.xl BWtest.F90 -W@,-v
#define SMEM_SIZE 32  
PROGRAM MAIN

      USE iso_c_binding
      TYPE compact_op1
         INTEGER(c_int) mo,no
         real(kind=c_double), dimension(:,:), allocatable :: art
         real(kind=c_double), dimension(:,:), allocatable :: ar
         real(kind=c_double), dimension(:,:), allocatable :: rc
      END TYPE compact_op1
      real(kind=c_double), allocatable, dimension(:,:,:) :: v
      real(kind=c_double), allocatable, dimension(:,:,:) :: dv
      real(kind=c_double), allocatable, dimension(:,:,:) :: vbr1, vbr2
      real(kind=c_double) :: scalefac
      real(kind=c_double) :: this_sum, this_v
      real(kind=c_double), allocatable, dimension(:,:) :: op_ar
      real(kind=c_double), dimension(-2:SMEM_SIZE+3) :: smem
      !real(kind=c_double), allocatable, dimension(-2:259) :: smem
      TYPE(compact_op1) :: op
      integer :: i,j,k,az,ay,ax,q,z

      allocate(v(262,262,262))
      allocate(dv(262,262,262))
      allocate(vbr1(3,262,262))
      allocate(vbr2(3,262,262))
      allocate(op_ar(7,262))
      allocate(op%ar(7,262))

      CALL srand(97035)
      CALL random_number(v)
      CALL random_number(vbr1)
      CALL random_number(vbr2)
      CALL random_number(op%ar)

      scalefac = 0.143*2
      ax=256
      ay=256
      az=256

      ! ****  from compact.f90
#ifndef CACHE
      do q=1,10
      !$acc parallel loop collapse(3) copyin(v,op,op%ar,vbr1,vbr2) copyout(dv)
#ifdef OLD
      !$omp target map(to:v,op%ar,vbr1,vbr2) map(from:dv)
      !$omp parallel do collapse(3) private(this_sum, this_v, z)
#else
      !$omp target teams distribute parallel do collapse(2) &
      !$omp map(to:v,op%ar,vbr1,vbr2,scalefac) map(from:dv) private(i,this_v,this_sum)
#endif
      do k=1,az
      do j=1,ay
      do i=1,ax
        this_sum = 0.0d0
        !$acc loop seq
        do z=-3,3
          if (z+i<1) then
             this_v = vbr1(z+i+3,j,k)
          elseif (z+i>ax) then
             this_v = vbr2(z+i-ax,j,k)
          else
             this_v = v(z+i,j,k)
          endif
          this_sum = this_sum + op%ar(z+4,i)*this_v
        end do
        dv(i,j,k) = this_sum * scalefac
      end do
      end do
      end do
#ifdef OLD
      !$omp end target
#endif
      v = dv
      end do

#else
      do q=1,10
      !$acc parallel loop gang collapse(3) copyin(v,op,op%ar,vbr1,vbr2) copyout(dv) private(smem) vector_length(SMEM_SIZE)
      !$omp target map(to:v,op%ar,vbr1,vbr2) map(from:dv)
      !!$omp teams distribute shared(smem) collapse(3) private(this_sum,this_v,z)
      !$omp parallel do collapse(3) private(this_sum, this_v, z) 
      do k=1,az
      do j=1,ay
      do i=0,ax,SMEM_SIZE
        !$acc cache(smem) 
        !$acc loop vector 
        do ii=1-3,SMEM_SIZE+3
          if (ii+i<1) then
             smem(ii) = vbr1(ii+i+3,j,k)
          elseif (ii+i>ax) then
             smem(ii) = vbr2(ii+i-ax,j,k)
          else
             smem(ii) = v(ii+i,j,k)
          endif
        end do
        !$acc loop vector
        do ii=1,SMEM_SIZE
        this_sum = 0.0d0
        !$acc loop seq
          do z=-3,3
            this_v = smem(z+ii)
            this_sum = this_sum + op%ar(z+4,ii+i)*this_v
          end do
        dv(ii+i,j,k) = this_sum * scalefac
        end do
      end do
      end do
      end do
      !$omp end target
      v = dv
      end do
#endif

      IF (abs(dv(104,37,214)  - 1.107543289) < 0.00001) THEN
         PRINT *, "PGI OK"
      ELSEIF (abs(dv(104,37,214) - 1.554295199) < 0.00001) THEN
         PRINT *, "XL OK"
      ELSE
         PRINT *, "********  ERROR  ******"
         PRINT *, v(104,37,214)
      ENDIF

END PROGRAM MAIN
