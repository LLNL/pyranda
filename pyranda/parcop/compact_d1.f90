module LES_compact_d1
  USE iso_c_binding
  USE MPI
  !USE LES_input, ONLY : bpp_lus_opt,use_ppent_opt,directcom,gpu_kernel
  USE LES_stencils
  USE LES_comm,  ONLY : mpierr,mpistatus,master
  use LES_pentadiagonal, ONLY :   ppentlus, &
    bpentLUS3x, ppentLUS3x, bpentLUS3y, ppentLUS3y, bpentLUS3z, ppentLUS3z, &
    btrid_block4_lus, ptrid_block4_lus, bpentLUS2y 
  use LES_compact_basetype, only : compact_op1, &
  																 vbr1x,vbr2x,vbs1x,vbs2x,vbr1y,vbr2y,vbs1y,vbs2y, &
  																 vbr1z,vbr2z,vbs1z,vbs2z,dvopx,dvox,dvopy,dvoy, &
  																 dvopz,dvoz

  IMPLICIT NONE

  integer :: directcom = 1
  integer :: gpu_kernel = 1
  logical(c_bool) :: bpp_lus_opt = .true. , use_ppent_opt = .true. 

  
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double, one=1.0_c_double
  LOGICAL(c_bool) :: debug=.false.

  ! these extensions override the generic operations
  
  type, extends(compact_op1) :: compact_op1_d1  ! custom operator (backward compat.)
  contains
    procedure :: evalx => eval_compact_op1x_d1  ! matrix version, OMP version with GPU support
    procedure :: evaly => eval_compact_op1y_d1  ! matrix version, OMP version with GPU support
    procedure :: evalz => eval_compact_op1z_d1  ! matrix version, OMP version with GPU support
  end type
 
contains

! d1 operations

  subroutine eval_compact_op1x_d1(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2)  ! generalized, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double), dimension(ay,ax) :: dv_tran
    integer :: nb,nsr,i,j,k
    integer :: nor,nir,nr,nol,nl,ni,np,dc ! surrogates
    character(len=160) :: filename
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
			do k=1,az
			 do j=1,ay
			  do i=1,ax
         dv(i,j,k) = zero
        end do ! k
       end do ! j
      end do ! i
      !$omp end target teams distribute parallel do
!      print *,'null op in x'
      RETURN
    END IF

    if( ax /= op%m ) then
      print *,'*** error: mismatch in x operation size ***',ax,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in x stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ; dc = op%directcom
    if( present(vb1) ) nb = size(vb1,1)

!---------------------------------------------------------------------------------------------------
! ghost data   
!---------------------------------------------------------------------------------------------------
    ! explicit part
    if( np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        vbs2x(1:nor,j,k) = v(ax-2:ax,j,k)
        vbs1x(1:nor,j,k) = v(1:3,j,k)
       end do ! j
      end do ! k
      !$omp end target teams distribute parallel do
      nsr = size(vbs1x)
      !$omp target data use_device_ptr(vbs1x,vbs2x,vbr1x,vbr2x) if(gpu_kernel==1)   ! use GPU Direct MPI
      call MPI_Sendrecv( vbs2x, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1x, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1x, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2x, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      !$omp end target data
    else if( op%periodic ) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        vbr1x(1:nor,j,k) = v(ax-2:ax,j,k)
        vbr2x(1:nor,j,k) = v(1:3,j,k)
       end do ! j
      end do ! k
      !$omp end target teams distribute parallel do
      if( debug ) print *,'periodic in x'
    endif
!---------------------------------------------------------------------------------------------------
! end ghost data transfer

!    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
!    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data

!---------------------------------------------------------------------------------------------------      
! Stencil Computations
    if( op%lo == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        dv(1,j,k) = sum(op%ar(4:7,1)*v(1:4,j,k))
        dv(2,j,k) = sum(op%ar(3:7,2)*v(1:5,j,k))
        dv(3,j,k) = sum(op%ar(2:7,3)*v(1:6,j,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        dv(1,j,k) = sum(op%ar(5:7,1)*(v(2:4,j,k)-v(1,j,k)))
        dv(2,j,k) = sum(op%ar(4:7,2)*(v(2:5,j,k)-v(1,j,k)))
        dv(3,j,k) = op%ar(5,3)*(v(4,j,k)-v(2,j,k))+op%ar(6,3)*(v(5,j,k)-v(1,j,k))+op%ar(7,3)*(v(6,j,k)-v(1,j,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
       dv(1,j,k) = op%ar(5,1)*(v(2,j,k)-vbr1x(3,j,k))+op%ar(6,1)*(v(3,j,k)-vbr1x(2,j,k))+op%ar(7,1)*(v(4,j,k)-vbr1x(1,j,k))
       dv(2,j,k) = op%ar(5,2)*(v(3,j,k)-v(1,j,k))+op%ar(6,2)*(v(4,j,k)-vbr1x(3,j,k))+op%ar(7,2)*(v(5,j,k)-vbr1x(2,j,k))
       dv(3,j,k) = op%ar(5,3)*(v(4,j,k)-v(2,j,k))+op%ar(6,3)*(v(5,j,k)-v(1,j,k))+op%ar(7,3)*(v(6,j,k)-vbr1x(3,j,k))
      end do
     end do
     !$omp end target teams distribute parallel do
    endif
    !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    do k=1,az
     do j=1,ay
      do i=4,ax-3
       dv(i,j,k) = op%ar(5,i)*(v(i+1,j,k)-v(i-1,j,k))+op%ar(6,i)*(v(i+2,j,k)-v(i-2,j,k))+op%ar(7,i)*(v(i+3,j,k)-v(i-3,j,k))
      end do
     end do
    end do
    !$omp end target teams distribute parallel do
    if( op%hi == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*v(ax-5:ax,j,k))
        dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*v(ax-4:ax,j,k))
        dv(ax,j,k)   = sum(op%ar(1:4,ax  )*v(ax-3:ax,j,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        dv(ax-2,j,k) = op%ar(1,ax-2)*(v(ax-5,j,k)-v(ax,j,k))+op%ar(2,ax-2)*(v(ax-4,j,k)-v(ax,j,k))+op%ar(3,ax-2)*(v(ax-3,j,k)-v(ax-1,j,k))
        dv(ax-1,j,k) = sum(op%ar(1:4,ax-1)*(v(ax-4:ax-1,j,k)-v(ax,j,k)))
        dv(ax,j,k)   = sum(op%ar(1:3,ax  )*(v(ax-3:ax-1,j,k)-v(ax,j,k)))
       end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(ax-2,j,k) = op%ar(5,ax-2)*(v(ax-1,j,k)-v(ax-3,j,k))+op%ar(6,ax-2)*(v(ax,j,k)-v(ax-4,j,k))+op%ar(7,ax-2)*(vbr2x(1,j,k)-v(ax-5,j,k))
       dv(ax-1,j,k) = op%ar(5,ax-1)*(v(ax,j,k)-v(ax-2,j,k))+op%ar(6,ax-1)*(vbr2x(1,j,k)-v(ax-3,j,k))+op%ar(7,ax-1)*(vbr2x(2,j,k)-v(ax-4,j,k))
       dv(ax,j,k)   = op%ar(5,ax  )*(vbr2x(1,j,k)-v(ax-1,j,k))+op%ar(6,ax)*(vbr2x(2,j,k)-v(ax-2,j,k))+op%ar(7,ax  )*(vbr2x(3,j,k)-v(ax-3,j,k))
      end do
     end do
     !$omp end target teams distribute parallel do
    endif
!---------------------------------------------------------------------------------------------------      
! Explicit Boundary Conditions
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then 
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(1,j,k)=dv1(j,k)   ! supply lower solution
      end do ! j
     end do ! k
     !$omp end target teams distribute parallel do
    end if
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then 
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(ax,j,k)=dv2(j,k)  ! supply upper solution
      end do ! j
     end do ! k
     !$omp end target teams distribute parallel do
    end if
    if( .not. op%implicit_op ) return
!---------------------------------------------------------------------------------------------------   
! implicit part
    if( np == 1 .and. op%periodic ) then			! No current GPU support for bpp_lus_opt
      if (bpp_lus_opt .and. (gpu_kernel==0) ) then
          do k = 1, az
             dv_tran = transpose(dv(:,:,k))
             call ppentlus_f77(ay,ax,1,op%al(1,1),dv_tran(1,1))
             dv(:,:,k) = transpose(dv_tran)
          end do
      else
         call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
      end if
    else
      if (bpp_lus_opt .and. (gpu_kernel==0) ) then
         do k = 1, az
            ! Often faster to transpose the data in 2D slabs and solve in y
            ! direction so that stride-1 X-direction will vectorize and expose
            ! more instruction level parallelism.
            dv_tran = transpose(dv(:,:,k))
            call bpentLUS2y(op%al,dv_tran,op%m,ay,ax)
            dv(:,:,k) = transpose(dv_tran)
         end do
      else
     		call bpentLUS3x(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
      endif ! bpp_lus_opt
    endif	! np 
    if( np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do j=1,ay
        dvopx(1:2,j,k) = dv(1:2,j,k)
        dvopx(3:4,j,k) = dv(ax-1:ax,j,k)
       end do ! k
      end do ! j
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then 
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az
        do j=1,ay
         dvopx(1:2,j,k) = zero
        end do ! k
       end do ! j
       !$omp end target teams distribute parallel do
      end if
      if( op%hi == MPI_PROC_NULL ) then 
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az
        do j=1,ay
         dvopx(3:4,j,k) = zero
        end do ! k
       end do ! j
       !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopx)
      select case( dc )
      case( 1 ) ! mpi_allgather
      	!$omp target data use_device_ptr(dvopx,dvox) if(gpu_kernel==1)
        call mpi_allgather(dvopx,nsr,MPI_DOUBLE_PRECISION,dvox,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        if( op%periodic ) then
         	call ptrid_block4_lus( op%aa, dvox, np, ay, az )	! cpu solver
        else
         	call btrid_block4_lus( op%aa, dvox, np, ay, az )
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
        	 !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvox(3,j,k,op%lo) - op%rc(i,2)*dvox(4,j,k,op%lo) &
                        &           - op%rc(i,3)*dvox(1,j,k,op%hi) - op%rc(i,4)*dvox(2,j,k,op%hi)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        else if (op%lo /= MPI_PROC_NULL) then
           !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvox(3,j,k,op%lo) - op%rc(i,2)*dvox(4,j,k,op%lo)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        else if (op%hi /= MPI_PROC_NULL) then
           !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,3)*dvox(1,j,k,op%hi) - op%rc(i,4)*dvox(2,j,k,op%hi)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        end if
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvopx,nsr,MPI_DOUBLE_PRECISION,dvox,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
           	call ptrid_block4_lus( op%aa, dvox, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvox, np, ay, az )
          endif
        else
          dvox = zero
        endif
        ! shuffle solution
        dvopx(3:4,:,:) = dvox(1:2,:,:,0)
        do i=0,np-2
          dvopx(1:2,:,:) = dvox(3:4,:,:,i)
          dvox(3:4,:,:,i) = dvox(1:2,:,:,i+1)
          dvox(1:2,:,:,i+1) = dvopx(1:2,:,:)
        end do
        dvox(1:2,:,:,0) = dvox(3:4,:,:,np-1)
        dvox(3:4,:,:,np-1) = dvopx(3:4,:,:)
        call mpi_scatter(dvox,nsr,MPI_DOUBLE_PRECISION,dvopx,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvopx(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvopx(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvopx(:,j,k))
      end select
    endif
!---------------------------------------------------------------------------------------------------      
  end subroutine eval_compact_op1x_d1

  subroutine eval_compact_op1y_d1(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    integer :: nb,nsr,i,j,k,n
    integer :: nor,nir,nr,nol,nl,ni,np,dc  ! surrogates
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
    	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
			do k=1,az
			 do j=1,ay
			  do i=1,ax    
         dv(i,j,k) = zero
        end do ! k
       end do ! j
      end do ! i
      !$omp end target teams distribute parallel do
      print *,'null op in y'
      RETURN
    END IF
    !ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)   ! Vestigial from CPU code
    if( ay /= op%m ) then
      print *,'*** error: mismatch in y operation size ***',ax,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in y stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr			! Set surrogates
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np  ; dc = op%directcom
    if( present(vb1) ) nb = size(vb1,2)
!---------------------------------------------------------------------------------------------------      
! Arrange halo data for MPI    
    ! explicit part
    if( op%np > 1 ) then  ! use parallel solver
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az
    		do i=1,ax
      		vbs2y(i,1:nor,k) = v(i,ay-2:ay,k)		! Fill vbs on GPU
      		vbs1y(i,1:nor,k) = v(i,1:3,k)
      	end do ! k
      end do ! i
      !$omp end target teams distribute parallel do
      nsr = size(vbs1y)
      !$omp target data use_device_ptr(vbs1y,vbs2y,vbr1y,vbr2y)	 if(gpu_kernel==1)			! Using Peer-to-peer memory transfers
      call MPI_Sendrecv( vbs2y, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1y, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1y, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2y, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      !$omp end target data
    else if( op%periodic ) then
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az
    	 do i=1,ax
      	vbr1y(i,1:nor,k) = v(i,ay-2:ay,k)
      	vbr2y(i,1:nor,k) = v(i,1:3,k)
       end do ! i
      end do ! k
      !$omp end target teams distribute parallel do
      if( debug ) print *,'periodic in y'
    endif
!---------------------------------------------------------------------------------------------------
! end ghost data transfer

!---------------------------------------------------------------------------------------------------      
! Stencil Computations
    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do i=1,ax
        dv(i,1,k) = sum(op%ar(4:7,1)*v(i,1:4,k))
        dv(i,2,k) = sum(op%ar(3:7,2)*v(i,1:5,k))
        dv(i,3,k) = sum(op%ar(2:7,3)*v(i,1:6,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do i=1,ax
        dv(i,1,k) = sum(op%ar(5:7,1)*(v(i,2:4,k)-v(i,1,k)))
        dv(i,2,k) = sum(op%ar(4:7,2)*(v(i,2:5,k)-v(i,1,k)))
        dv(i,3,k) = op%ar(5,3)*(v(i,4,k)-v(i,2,k))+op%ar(6,3)*(v(i,5,k)-v(i,1,k))+op%ar(7,3)*(v(i,6,k)-v(i,1,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do i=1,ax
       dv(i,1,k) = op%ar(5,1)*(v(i,2,k)-vbr1y(i,3,k))+op%ar(6,1)*(v(i,3,k)-vbr1y(i,2,k))+op%ar(7,1)*(v(i,4,k)-vbr1y(i,1,k))
       dv(i,2,k) = op%ar(5,2)*(v(i,3,k)-v(i,1,k))+op%ar(6,2)*(v(i,4,k)-vbr1y(i,3,k))+op%ar(7,2)*(v(i,5,k)-vbr1y(i,2,k))
       dv(i,3,k) = op%ar(5,3)*(v(i,4,k)-v(i,2,k))+op%ar(6,3)*(v(i,5,k)-v(i,1,k))+op%ar(7,3)*(v(i,6,k)-vbr1y(i,3,k))
      end do
      end do
      !$omp end target teams distribute parallel do
    endif
     !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
     do k=1,az
      do j=4,ay-3
       do i=1,ax
        dv(i,j,k) = op%ar(5,j)*(v(i,j+1,k)-v(i,j-1,k))+op%ar(6,j)*(v(i,j+2,k)-v(i,j-2,k))+op%ar(7,j)*(v(i,j+3,k)-v(i,j-3,k))
       end do
      end do
     end do
     !$omp end target teams distribute parallel do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do i=1,ax
        dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*v(i,ay-5:ay,k))
        dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*v(i,ay-4:ay,k))
        dv(i,ay,k)   = sum(op%ar(1:4,ay  )*v(i,ay-3:ay,k))
       end do
      end do
      !$omp end target teams distribute parallel do
     else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do i=1,ax
        dv(i,ay-2,k) = op%ar(1,ay-2)*(v(i,ay-5,k)-v(i,ay,k))+op%ar(2,ay-2)*(v(i,ay-4,k)-v(i,ay,k))+op%ar(3,ay-2)*(v(i,ay-3,k)-v(i,ay-1,k))
        dv(i,ay-1,k) = sum(op%ar(1:4,ay-1)*(v(i,ay-4:ay-1,k)-v(i,ay,k)))
        dv(i,ay,k)   = sum(op%ar(1:3,ay  )*(v(i,ay-3:ay-1,k)-v(i,ay,k)))
       end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do i=1,ax
       dv(i,ay-2,k) = op%ar(5,ay-2)*(v(i,ay-1,k)-v(i,ay-3,k))+op%ar(6,ay-2)*(v(i,ay,k)-v(i,ay-4,k))+op%ar(7,ay-2)*(vbr2y(i,1,k)-v(i,ay-5,k))
       dv(i,ay-1,k) = op%ar(5,ay-1)*(v(i,ay,k)-v(i,ay-2,k))+op%ar(6,ay-1)*(vbr2y(i,1,k)-v(i,ay-3,k))+op%ar(7,ay-1)*(vbr2y(i,2,k)-v(i,ay-4,k))
       dv(i,ay,k)   = op%ar(5,ay  )*(vbr2y(i,1,k)-v(i,ay-1,k))+op%ar(6,ay)*(vbr2y(i,2,k)-v(i,ay-2,k))+op%ar(7,ay  )*(vbr2y(i,3,k)-v(i,ay-3,k))
      end do
     end do
     !$omp end target teams distribute parallel do
    endif
!---------------------------------------------------------------------------------------------------      
! Explicit Boundary Conditions
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do i=1,ax
       dv(i,1,k)=dv1(i,k)     ! supply lower solution
      end do ! i
     end do ! k
     !$omp end target teams distribute parallel do
    end if
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do i=1,ax
       dv(i,ay+1,k)=dv2(i,k)  ! supply upper solution
			end do !i
     end do !k
     !$omp end target teams distribute parallel do
    end if
    if( .not. op%implicit_op ) return
!---------------------------------------------------------------------------------------------------      
! Solve implicit portion - pentadiagonal solve
    if( op%np == 1 .and. op%periodic ) then
     if ( bpp_lus_opt .and. (gpu_kernel==0)) then
      call ppentlus_f77(ax,ay,az,op%al(1,1),dv(1,1,1))  ! periodic solution on a single process
     else
      call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
     end if
    else
     call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif	! np
    if( op%np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
       do i=1,ax
        dvopy(1:2,i,k) = dv(i,1:2,k)
        dvopy(3:4,i,k) = dv(i,ay-1:ay,k)
       end do
      end do
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az
        do i=1,ax
         dvopy(1:2,i,k) = zero
        end do ! i
       end do ! k
       !$omp end target teams distribute parallel do
      end if
      if( op%hi == MPI_PROC_NULL ) then 
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az
        do i=1,ax
         dvopy(3:4,i,k) = zero
        end do ! i
       end do ! k
       !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopy)
!---------------------------------------------------------------------------------------------------      
! Global Solve
      select case( dc )
      case( 1 ) ! mpi_allgather
      	!$omp target data use_device_ptr(dvopy, dvoy) if(gpu_kernel==1)
        call mpi_allgather(dvopy,nsr,MPI_DOUBLE_PRECISION,dvoy,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        
        if( op%periodic ) then
         	call ptrid_block4_lus( op%aa, dvoy, op%np, ax, az )
        else
         	call btrid_block4_lus( op%aa, dvoy, op%np, ax, az ) ! CPU solver
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
         !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
         do k = 1, az
          do j = 1, ay
           do i = 1, ax
            dv(i,j,k) = dv(i,j,k) - op%rc(j,1)*dvoy(3,i,k,op%lo) - op%rc(j,2)*dvoy(4,i,k,op%lo) &
                        &           - op%rc(j,3)*dvoy(1,i,k,op%hi) - op%rc(j,4)*dvoy(2,i,k,op%hi)
           end do
          end do
         end do
         !$omp end target teams distribute parallel do
        else if (op%lo /= MPI_PROC_NULL) then
         !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
         do k = 1, az
          do j = 1, ay
           do i = 1, ax
            dv(i,j,k) = dv(i,j,k) - op%rc(j,1)*dvoy(3,i,k,op%lo) - op%rc(j,2)*dvoy(4,i,k,op%lo)
           end do
          end do
         end do
         !$omp end target teams distribute parallel do
        else if (op%hi /= MPI_PROC_NULL) then
         !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
         do k = 1, az
          do j = 1, ay
           do i = 1, ax
            dv(i,j,k) = dv(i,j,k) - op%rc(j,3)*dvoy(1,i,k,op%hi) - op%rc(j,4)*dvoy(2,i,k,op%hi)
           end do
          end do
         end do
         !$omp end target teams distribute parallel do
        end if
        
      case( 2 ) ! mpi_gather/scatter
				!$omp target data use_device_ptr(dvopy, dvoy) if(gpu_kernel==1)
        call mpi_gather(dvopy,nsr,MPI_DOUBLE_PRECISION,dvoy,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        !$omp end target data
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvoy, op%np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvoy, op%np, ax, az )
          endif
        else
         !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
         do k=1,az
        	do i=1,ax
           dvoy(:,i,k,:) = zero
          end do ! i
         end do ! k
         !$omp end target teams distribute parallel do
        endif
        ! shuffle solution
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do k=1,az
         do i=1,ax
          dvopy(3:4,i,k) = dvoy(1:2,i,k,0)
         end do ! i
        end do ! k
        !$omp end target teams distribute parallel do
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do k=1,az
         do i=1,ax
        	do n=0,op%np-2
           dvopy(1:2,i,k) = dvoy(3:4,i,k,n)
           dvoy(3:4,i,k,n) = dvoy(1:2,i,k,n+1)
           dvoy(1:2,i,k,n+1) = dvopy(1:2,i,k)
          end do ! np
         end do ! i
        end do ! k
        !$omp end target teams distribute parallel do
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do k=1,az
         do i=1,ax
          dvoy(1:2,i,k,0) = dvoy(3:4,i,k,op%np-1)
          dvoy(3:4,i,k,op%np-1) = dvopy(3:4,i,k)
         end do ! i
        end do ! k
        !$omp end target teams distribute parallel do

        !$omp target data use_device_ptr(dvoy,dvopy) if(gpu_kernel==1)
        call mpi_scatter(dvoy,nsr,MPI_DOUBLE_PRECISION,dvopy,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        !$omp end target data
        if( op%lo == MPI_PROC_NULL ) then
         !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	   do k=1,az
    	    do i=1,ax
    	     dvopy(1:2,i,k) = zero
    	    end do ! i
         end do ! k
         !$omp end target teams distribute parallel do
        end if
        if( op%hi == MPI_PROC_NULL ) then
         !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	   do k=1,az
    	    do i=1,ax
           dvopy(3:4,i,k) = zero
          end do ! i
         end do ! k
         !$omp end target teams distribute parallel do
        end if
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	   do k=1,az
    	    do j=1,ay
    	     do i=1,ax
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvopy(1:ni,i,k))
           end do ! i
          end do ! j
         end do ! k
         !$omp end target teams distribute parallel do
      end select
    endif ! np
!---------------------------------------------------------------------------------------------------      
  end subroutine eval_compact_op1y_d1
    
  subroutine eval_compact_op1z_d1(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2)  ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    integer :: nb,nsr,i,j,k,n
    integer :: nor,nir,nr,nol,nl,ni,np,dc  ! surrogates
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
     !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
			do k=1,az
			 do j=1,ay
			  do i=1,ax  
         dv(i,j,k) = zero
        end do ! i
       end do ! j
      end do ! k
      !$omp end target teams distribute parallel do
!      print *,'null op in z'
      RETURN
    END IF
    if( az /= op%m ) then
      print *,'*** error: mismatch in z operation size ***',ax,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in z stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ; dc=op%directcom
    if( present(vb1) ) nb = size(vb1,3)
!---------------------------------------------------------------------------------------------------      
! Arrange halo data for MPI    
    ! explicit part
    if( np > 1 ) then  ! use parallel solver
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) ! Fill vbs on GPU
		  do j=1,ay
		   do i=1,ax
		    vbs2z(i,j,1:nor) = v(i,j,az-2:az)
		    vbs1z(i,j,1:nor) = v(i,j,1:3)
		   end do ! i
		  end do ! j
		  !$omp end target teams distribute parallel do
      nsr = size(vbs1z)
      !$omp target data use_device_ptr(vbs1z,vbs2z,vbr1z,vbr2z) if(gpu_kernel==1)
      call MPI_Sendrecv( vbs2z, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1z, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1z, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2z, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      !$omp end target data
    else if( op%periodic ) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
      do i=1,ax
       vbr1z(i,j,1:nor) = v(i,j,az-2:az)
       vbr2z(i,j,1:nor) = v(i,j,1:3)
      end do ! i
     end do ! j
     !$omp end target teams distribute parallel do
      if( debug ) print *,'periodic in z'
    endif
!---------------------------------------------------------------------------------------------------      
!---------------------------------------------------------------------------------------------------
! end ghost data transfer 

    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = sum(op%ar(4:7,1)*v(i,j,1:4))
      dv(i,j,2) = sum(op%ar(3:7,2)*v(i,j,1:5))
      dv(i,j,3) = sum(op%ar(2:7,3)*v(i,j,1:6))
     end do
     end do
     !$omp end target teams distribute parallel do
     else
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = sum(op%ar(5:7,1)*(v(i,j,2:4)-v(i,j,1)))
      dv(i,j,2) = sum(op%ar(4:7,2)*(v(i,j,2:5)-v(i,j,1)))
      dv(i,j,3) = op%ar(5,3)*(v(i,j,4)-v(i,j,2))+op%ar(6,3)*(v(i,j,5)-v(i,j,1))+op%ar(7,3)*(v(i,j,6)-v(i,j,1))
     end do
     end do
     !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = op%ar(5,1)*(v(i,j,2)-vbr1z(i,j,3))+op%ar(6,1)*(v(i,j,3)-vbr1z(i,j,2))+op%ar(7,1)*(v(i,j,4)-vbr1z(i,j,1))
      dv(i,j,2) = op%ar(5,2)*(v(i,j,3)-v(i,j,1))+op%ar(6,2)*(v(i,j,4)-vbr1z(i,j,3))+op%ar(7,2)*(v(i,j,5)-vbr1z(i,j,2))
      dv(i,j,3) = op%ar(5,3)*(v(i,j,4)-v(i,j,2))+op%ar(6,3)*(v(i,j,5)-v(i,j,1))+op%ar(7,3)*(v(i,j,6)-vbr1z(i,j,3))
     end do
     end do
     !$omp end target teams distribute parallel do
    endif
    !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    do k=4,az-3
    do j=1,ay
    do i=1,ax
      dv(i,j,k) = op%ar(5,k)*(v(i,j,k+1)-v(i,j,k-1))+op%ar(6,k)*(v(i,j,k+2)-v(i,j,k-2))+op%ar(7,k)*(v(i,j,k+3)-v(i,j,k-3))
    end do
    end do
    end do
    !$omp end target teams distribute parallel do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = sum(op%ar(1:6,az-2)*v(i,j,az-5:az))
      dv(i,j,az-1) = sum(op%ar(1:5,az-1)*v(i,j,az-4:az))
      dv(i,j,az)   = sum(op%ar(1:4,az  )*v(i,j,az-3:az))
     end do
     end do
     !$omp end target teams distribute parallel do
     else
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = op%ar(1,az-2)*(v(i,j,az-5)-v(i,j,az))+op%ar(2,az-2)*(v(i,j,az-4)-v(i,j,az))+op%ar(3,az-2)*(v(i,j,az-3)-v(i,j,az-1))
      dv(i,j,az-1) = sum(op%ar(1:4,az-1)*(v(i,j,az-4:az-1)-v(i,j,az)))
      dv(i,j,az)   = sum(op%ar(1:3,az  )*(v(i,j,az-3:az-1)-v(i,j,az)))
     end do
     end do
     !$omp end target teams distribute parallel do
     endif
    else ! centered interior weights
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = op%ar(5,az-2)*(v(i,j,az-1)-v(i,j,az-3))+op%ar(6,az-2)*(v(i,j,az)-v(i,j,az-4))+op%ar(7,az-2)*(vbr2z(i,j,1)-v(i,j,az-5))
      dv(i,j,az-1) = op%ar(5,az-1)*(v(i,j,az)-v(i,j,az-2))+op%ar(6,az-1)*(vbr2z(i,j,1)-v(i,j,az-3))+op%ar(7,az-1)*(vbr2z(i,j,2)-v(i,j,az-4))
      dv(i,j,az)   = op%ar(5,az  )*(vbr2z(i,j,1)-v(i,j,az-1))+op%ar(6,az)*(vbr2z(i,j,2)-v(i,j,az-2))+op%ar(7,az  )*(vbr2z(i,j,3)-v(i,j,az-3))
     end do
     end do
     !$omp end target teams distribute parallel do
    endif
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
      do i=1,ax
       dv(i,j,1)=dv1(i,j)   ! supply lower solution
      end do ! i
     end do ! j
     !$omp end target teams distribute parallel do
    end if
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do j=1,ay
      do i=1,ax
       dv(i,j,az)=dv2(i,j)  ! supply upper solution
      end do ! i
     end do ! j
     !$omp end target teams distribute parallel do
    end if
    if( .not. op%implicit_op ) return
    ! implicit part
    if( op%np == 1 .and. op%periodic ) then
       if (bpp_lus_opt .and. (gpu_kernel==0) ) then
          call ppentlus(3,use_ppent_opt,op%al,dv)
       else
          call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
       end if
    else
     	call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif	! np

    if( op%np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
       do i=1,ax
      	do k=1,2
         dvopz(k,i,j) = dv(i,j,k)
         dvopz(2+k,i,j) = dv(i,j,az-2+k)
      	end do ! k
       end do ! i
      end do ! j
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then
		    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
		    do j=1,ay
		     do i=1,ax
		    	dvopz(1:2,i,j) = zero
		     end do ! j
		    end do ! i
		    !$omp end target teams distribute parallel do
		  end if
		  if( op%hi == MPI_PROC_NULL ) then 
		    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
		    do j=1,ay
		     do i=1,ax
		      dvopz(3:4,i,j) = zero
		     end do ! i
		    end do ! j
		    !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopz)
      select case( dc )
      case( 1 ) ! mpi_allgather
      	!$omp target data use_device_ptr(dvopz, dvoz) if(gpu_kernel==1)
        call mpi_allgather(dvopz,nsr,MPI_DOUBLE_PRECISION,dvoz,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        if( op%periodic ) then
         	call ptrid_block4_lus( op%aa, dvoz, op%np, ax, ay ) ! cpu solver
        else
         	call btrid_block4_lus( op%aa, dvoz, op%np, ax, ay ) ! CPU solve
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
        	 !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,1)*dvoz(3,i,j,op%lo) - op%rc(k,2)*dvoz(4,i,j,op%lo) &
                        &           - op%rc(k,3)*dvoz(1,i,j,op%hi) - op%rc(k,4)*dvoz(2,i,j,op%hi)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        else if (op%lo /= MPI_PROC_NULL) then
        	 !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,1)*dvoz(3,i,j,op%lo) - op%rc(k,2)*dvoz(4,i,j,op%lo)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        else if (op%hi /= MPI_PROC_NULL) then
        	 !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,3)*dvoz(1,i,j,op%hi) - op%rc(k,4)*dvoz(2,i,j,op%hi)
           end do
           end do
           end do
           !$omp end target teams distribute parallel do
        end if
      case( 2 ) ! mpi_gather/scatter
      	!$omp target data use_device_ptr( dvopz,dvoz ) if(gpu_kernel==1)
        call mpi_gather(dvopz,nsr,MPI_DOUBLE_PRECISION,dvoz,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        !$omp end target data
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvoz, op%np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvoz, op%np, ax, ay )
          endif
        else
         !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
         do j=1,ay
          do i=1,ax
           dvoz(:,i,j,:) = zero
          end do ! j
         end do ! i
         !$omp end target teams distribute parallel do
        endif
        ! shuffle solution
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do j=1,ay
         do i=1,ax
          dvopz(3:4,i,j) = dvoz(1:2,i,j,0)
          do n=0,op%np-2
           dvopz(1:2,i,j) = dvoz(3:4,i,j,n)
           dvoz(3:4,i,j,n) = dvoz(1:2,i,j,n+1)
           dvoz(1:2,i,j,n+1) = dvopz(1:2,i,j)
          end do
          dvoz(1:2,i,j,0) = dvoz(3:4,i,j,op%np-1)
          dvoz(3:4,i,j,op%np-1) = dvopz(3:4,i,j)
         end do ! j
        end do ! i
        !$omp end target teams distribute parallel do
        !$omp target data use_device_ptr(dvoz,dvopz) if(gpu_kernel==1)
        call mpi_scatter(dvoz,nsr,MPI_DOUBLE_PRECISION,dvopz,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        !$omp end target data
        !$omp target teams if(gpu_kernel==1)
        if( op%lo == MPI_PROC_NULL ) then
         !$omp distribute parallel do collapse(2)
         do j=1,ay
          do i=1,ax
           dvopz(1:2,i,j) = zero
          end do ! j
         end do ! i
         !$omp end distribute parallel do
        end if
        if( op%hi == MPI_PROC_NULL ) then
        !$omp distribute parallel do collapse(2)
         do j=1,ay
          do i=1,ax
           dvopz(3:4,:,:) = zero
          end do ! j
         end do ! i
         !$omp end distribute parallel do
        end if
        !$omp distribute parallel do collapse(3)
         do k=1,az
          do j=1,ay
           do i=1,ax
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvopz(1:ni,i,j))
           end do ! i
          end do ! j
         end do ! k
         !$omp end distribute parallel do
         !$omp end target teams
      end select
      
    endif
  end subroutine eval_compact_op1z_d1


end module LES_compact_d1

