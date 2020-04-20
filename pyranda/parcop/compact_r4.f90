module LES_compact_r4
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
  
  type, extends(compact_op1) :: compact_op1_r4  ! custom operator (backward compat.)
  contains
    procedure :: evalx => eval_compact_op1x_r4  ! rhs stencil = 4, OMP version with GPU support
    procedure :: evaly => eval_compact_op1y_r4  ! rhs stencil = 4, OMP version with GPU support
    procedure :: evalz => eval_compact_op1z_r4  ! rhs stencil = 4, OMP version with GPU support
  end type

contains

! "optimized" r4 operators for backward compatibility with matrix.f

  subroutine eval_compact_op1x_r4(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2)  ! generalized, uses 1D op type
    implicit none
    class(compact_op1_r4), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double) :: vc
    integer :: nb,nsr,i,j,k,l
    integer :: nor,nir,nr,nol,nl,ni,np,dc  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) then
      	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
      	do k=1,az; do j=1,ay; do i=1,ax
      		dv(i,j,k) = zero
      	end do; end do; end do
      	!$omp end target teams distribute parallel do
      end if
      
      if( op%null_option == 1 ) then
      	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
      	do k=1,az; do j=1,ay; do i=1,ax
      		dv(i,j,k) = v(i,j,k)  ! filter / interp
				end do; end do; end do
				!$omp end target teams distribute parallel do
			end if
!      print *,'null op in x'
      RETURN
    END IF
    
    if( ax /= op%m ) then
      print *,'*** error: mismatch in x operation size ***',ax,op%m
      stop
    endif
    if( op%nor /= 4 ) then
      print *,'*** error: mismatch in x stencil size ***',4,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ; dc = op%directcom
    ! explicit part
    ! ghost data

    if( np > 1 ) then  ! use parallel solver
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az; do j=1,ay
      	vbs2x(1:nor,j,k) = v(ax-3:ax,j,k)
      	vbs1x(1:nor,j,k) = v(1:4,j,k)
      end do; end do
      !$omp end target teams distribute parallel do
      nsr = size(vbs1x)
      !$omp target data use_device_ptr(vbs1x,vbs2x,vbr1x,vbr2x) if(gpu_kernel==1)
      call MPI_Sendrecv( vbs2x, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1x, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1x, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2x, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      !$omp end target data
    else if( op%periodic ) then
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az; do j=1,ay
      	vbr1x(1:nor,j,k) = v(ax-3:ax,j,k)
      	vbr2x(1:nor,j,k) = v(1:4,j,k)
      end do; end do
      !$omp end target teams distribute parallel do
    endif
        
    if( op%lo == MPI_PROC_NULL ) then
      if( op%bc(1) == -1 ) then ! antisymmetric BC, sumr != 0
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az
    	 do j=1,ay
          dv(1,j,k) = sum(op%ar(5:9,1)*v(1:5,j,k))
          dv(2,j,k) = sum(op%ar(4:9,2)*v(1:6,j,k))
          dv(3,j,k) = sum(op%ar(3:9,3)*v(1:7,j,k))
          dv(4,j,k) = sum(op%ar(2:9,4)*v(1:8,j,k))
    	 end do
    	end do
    	!$omp end target teams distribute parallel do
      else ! symmetric/one-sided BC, assume sumr=0
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az
    	 do j=1,ay
          dv(1,j,k) = sum(op%ar(6:9,1)*(v(2:5,j,k)-v(1,j,k)))
          dv(2,j,k) = sum(op%ar(5:9,2)*(v(2:6,j,k)-v(1,j,k)))
          dv(3,j,k) = sum(op%ar(4:9,3)*(v(2:7,j,k)-v(1,j,k)))
          dv(4,j,k) = sum(op%ar(3:9,4)*(v(2:8,j,k)-v(1,j,k)))
    	 end do
    	end do
    	!$omp end target teams distribute parallel do
      endif
    else
    	!$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    	do k=1,az
    	 do j=1,ay
          dv(1,j,k) = sum(op%ar(1:4,1)*(vbr1x(1:4,j,k)-v(1,j,k)))+sum(op%ar(5:9,1)*(v(1:5,j,k)-v(1,j,k)))
          dv(2,j,k) = sum(op%ar(1:3,2)*(vbr1x(2:4,j,k)-v(2,j,k)))+sum(op%ar(4:9,2)*(v(1:6,j,k)-v(2,j,k)))
          dv(3,j,k) = sum(op%ar(1:2,3)*(vbr1x(3:4,j,k)-v(3,j,k)))+sum(op%ar(3:9,3)*(v(1:7,j,k)-v(3,j,k)))
          dv(4,j,k) = sum(op%ar(1:1,4)*(vbr1x(4:4,j,k)-v(4,j,k)))+sum(op%ar(2:9,4)*(v(1:8,j,k)-v(4,j,k)))
    	 end do
    	end do
    	!$omp end target teams distribute parallel do
    endif

    !$omp target teams distribute parallel do collapse(3) private(l,vc) if(gpu_kernel==1)
    do k=1,az
      do j=1,ay
        do i=5,ax-4
          dv(i,j,k) = zero
          vc = v(i,j,k)
          do l=1,9
             dv(i,j,k) = dv(i,j,k) + ( v(l+i-5,j,k) - vc )*op%ar(l,i)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! antisymmetric BC, sumr != 0
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(ax-3,j,k) = sum(op%ar(1:8,ax-3)*v(ax-7:ax,j,k))
       dv(ax-2,j,k) = sum(op%ar(1:7,ax-2)*v(ax-6:ax,j,k))
       dv(ax-1,j,k) = sum(op%ar(1:6,ax-1)*v(ax-5:ax,j,k))
       dv(ax  ,j,k) = sum(op%ar(1:5,ax  )*v(ax-4:ax,j,k))
      end do
     end do
     !$omp end target teams distribute parallel do
     else ! symmetric/one-sided BC, assume sumr=0
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(ax-3,j,k) = sum(op%ar(1:7,ax-3)*(v(ax-7:ax-1,j,k)-v(ax,j,k)))
       dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*(v(ax-6:ax-1,j,k)-v(ax,j,k)))
       dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*(v(ax-5:ax-1,j,k)-v(ax,j,k)))
       dv(ax  ,j,k) = sum(op%ar(1:4,ax  )*(v(ax-4:ax-1,j,k)-v(ax,j,k)))
      end do
     end do
     !$omp end target teams distribute parallel do
     endif
    else
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az
      do j=1,ay
       dv(ax-3,j,k) = sum(op%ar(1:8,ax-3)*(v(ax-7:ax,j,k)-v(ax-3,j,k)))+sum(op%ar(9:9,ax-3)*(vbr2x(1:1,j,k)-v(ax-3,j,k)))
       dv(ax-2,j,k) = sum(op%ar(1:7,ax-2)*(v(ax-6:ax,j,k)-v(ax-2,j,k)))+sum(op%ar(8:9,ax-2)*(vbr2x(1:2,j,k)-v(ax-2,j,k)))
       dv(ax-1,j,k) = sum(op%ar(1:6,ax-1)*(v(ax-5:ax,j,k)-v(ax-1,j,k)))+sum(op%ar(7:9,ax-1)*(vbr2x(1:3,j,k)-v(ax-1,j,k)))
       dv(ax  ,j,k) = sum(op%ar(1:5,ax  )*(v(ax-4:ax,j,k)-v(ax  ,j,k)))+sum(op%ar(6:9,ax  )*(vbr2x(1:4,j,k)-v(ax  ,j,k)))
      end do
     end do
     !$omp end target teams distribute parallel do
    endif
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az; do j=1,ay
      dv(1,j,k)=dv1(j,k)   ! supply lower solution
     end do; end do
     !$omp end target teams distribute parallel do
    endif
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then 
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
     do k=1,az; do j=1,ay
      dv(ax,j,k)=dv2(j,k)  ! supply upper solution
     end do; end do
     !$omp end target teams distribute parallel do
    endif
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) then 
       !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
       do k=1,az; do j=1,ay; do i=1,ax
        dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
       end do; end do; end do
       !$omp end target teams distribute parallel do
      end if
      return
    endif
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) then 
       !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
       do k=1,az; do j=1,ay; do i=1,ax
        dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
       end do; end do; end do
       !$omp end target teams distribute parallel do
      end if
      return
    endif
    ! parallel solver
    if( np > 1 ) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do j=1,ay;
       do i=1,2
        dvopx(i,j,k) = dv(i,j,k)
        dvopx(i+2,j,k) = dv(ax+i-2,j,k)
       end do;
      end do; end do
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az; do j=1,ay
        dvopx(1:2,j,k) = zero
       end do; end do
       !$omp end target teams distribute parallel do
      end if
      if( op%hi == MPI_PROC_NULL ) then 
       !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
       do k=1,az; do j=1,ay
        dvopx(3:4,j,k) = zero
       end do; end do
       !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopx)
      select case( dc )
      case( 1 ) ! mpi_allgather
        !$omp target data use_device_ptr(dvopx,dvox) if(gpu_kernel==1)
        call mpi_allgather(dvopx,nsr,MPI_DOUBLE_PRECISION,dvox,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        if( op%periodic ) then
         	call ptrid_block4_lus( op%aa, dvox, np, ay, az )	! CPU solve
        else
         	call btrid_block4_lus( op%aa, dvox, np, ay, az )	! CPU solve
        endif
        if( op%lo /= MPI_PROC_NULL ) then
         !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
         do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:2)*dvox(3:4,j,k,op%lo))
         end do; end do; end do
         !$omp end target teams distribute parallel do
        end if
        if( op%hi /= MPI_PROC_NULL ) then
         !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
         do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,3:4)*dvox(1:2,j,k,op%hi))
         end do; end do; end do
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
      
      if( op%null_option == 1 ) then
       !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
       do k=1,az; do j=1,ay; do i=1,ax
        dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
       end do; end do; end do
       !$omp end target teams distribute parallel do
      end if
    end if
  end subroutine eval_compact_op1x_r4


  subroutine eval_compact_op1y_r4(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2) ! nor=4, nol=2, uses 1D op type
    implicit none
    integer, intent(in) :: ax,ay,az
    class(compact_op1_r4), intent(in) :: op
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double) :: vc
    integer :: nb,nsr,i,j,k,l
    integer :: nor,nir,nr,nol,nl,ni,np,dc  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k) = zero
        end do; end do; end do
        !$omp end target teams distribute parallel do
      end if
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k) = v(i,j,k)  ! filter / interp
        end do; end do; end do
			  !$omp end target teams distribute parallel do
			end if
!      print *,'null op in y'
      RETURN
    END IF
    if( ay /= op%m ) then
      print *,'*** error: mismatch in y operation size ***',ay,op%m
      stop
    endif
    if( op%nor /= 4 ) then
      print *,'*** error: mismatch in y stencil size ***',4,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np  ; dc = op%directcom
    ! explicit part
    ! ghost data
    if( np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do i=1,ax
        vbs2y(i,1:nor,k) = v(i,ay-3:ay,k)
        vbs1y(i,1:nor,k) = v(i,1:4,k)
      end do; end do;
      !$omp end target teams distribute parallel do
      nsr = size(vbs1y)
      !$omp target data use_device_ptr(vbs1y,vbs2y,vbr1y,vbr2y) if(gpu_kernel==1)
      call MPI_Sendrecv( vbs2y, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1y, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1y, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2y, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      !$omp end target data
    else if( op%periodic ) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do i=1,ax
        vbr1y(i,1:nor,k) = v(i,ay-3:ay,k)
        vbr2y(i,1:nor,k) = v(i,1:4,k)
      end do; end do;
      !$omp end target teams distribute parallel do
      if( debug ) print *,'periodic in y'
    endif

    
    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! antisymmetric BC, sumr != 0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,1,k) = sum(op%ar(5:9,1)*v(i,1:5,k))
          dv(i,2,k) = sum(op%ar(4:9,2)*v(i,1:6,k))
          dv(i,3,k) = sum(op%ar(3:9,3)*v(i,1:7,k))
          dv(i,4,k) = sum(op%ar(2:9,4)*v(i,1:8,k))
        end do
      end do
      !$omp end target teams distribute parallel do
     else ! symmetric/one-sided BC, assume sumr=0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,1,k) = sum(op%ar(6:9,1)*(v(i,2:5,k)-v(i,1,k)))
          dv(i,2,k) = sum(op%ar(5:9,2)*(v(i,2:6,k)-v(i,1,k)))
          dv(i,3,k) = sum(op%ar(4:9,3)*(v(i,2:7,k)-v(i,1,k)))
          dv(i,4,k) = sum(op%ar(3:9,4)*(v(i,2:8,k)-v(i,1,k)))
        end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,1,k) = sum(op%ar(1:4,1)*(vbr1y(i,1:4,k)-v(i,1,k)))+sum(op%ar(5:9,1)*(v(i,1:5,k)-v(i,1,k)))
          dv(i,2,k) = sum(op%ar(1:3,2)*(vbr1y(i,2:4,k)-v(i,2,k)))+sum(op%ar(4:9,2)*(v(i,1:6,k)-v(i,2,k)))
          dv(i,3,k) = sum(op%ar(1:2,3)*(vbr1y(i,3:4,k)-v(i,3,k)))+sum(op%ar(3:9,3)*(v(i,1:7,k)-v(i,3,k)))
          dv(i,4,k) = sum(op%ar(1:1,4)*(vbr1y(i,4:4,k)-v(i,4,k)))+sum(op%ar(2:9,4)*(v(i,1:8,k)-v(i,4,k)))
        end do
      end do
      !$omp end target teams distribute parallel do
    endif
    !$omp target teams distribute parallel do collapse(3) private(l,vc) if(gpu_kernel==1)
    do k=1,az
      do j=5,ay-4
        do i=1,ax
          dv(i,j,k) = zero
          vc = v(i,j,k)
          do l=1,9
            dv(i,j,k) = dv(i,j,k) + ( v(i,l+j-5,k) - vc )*op%ar(l,j)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! antisymmetric BC, sumr != 0
     !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,ay-3,k) = sum(op%ar(1:8,ay-3)*v(i,ay-7:ay,k))
          dv(i,ay-2,k) = sum(op%ar(1:7,ay-2)*v(i,ay-6:ay,k))
          dv(i,ay-1,k) = sum(op%ar(1:6,ay-1)*v(i,ay-5:ay,k))
          dv(i,ay  ,k) = sum(op%ar(1:5,ay  )*v(i,ay-4:ay,k))
        end do
      end do
      !$omp end target teams distribute parallel do
     else ! symmetric/one-sided BC, assume sumr=0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,ay-3,k) = sum(op%ar(1:7,ay-3)*(v(i,ay-7:ay-1,k)-v(i,ay,k)))
          dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*(v(i,ay-6:ay-1,k)-v(i,ay,k)))
          dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*(v(i,ay-5:ay-1,k)-v(i,ay,k)))
          dv(i,ay  ,k) = sum(op%ar(1:4,ay  )*(v(i,ay-4:ay-1,k)-v(i,ay,k)))
        end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az
        do i=1,ax
          dv(i,ay-3,k) = sum(op%ar(1:8,ay-3)*(v(i,ay-7:ay,k)-v(i,ay-3,k)))+sum(op%ar(9:9,ay-3)*(vbr2y(i,1:1,k)-v(i,ay-3,k)))
          dv(i,ay-2,k) = sum(op%ar(1:7,ay-2)*(v(i,ay-6:ay,k)-v(i,ay-2,k)))+sum(op%ar(8:9,ay-2)*(vbr2y(i,1:2,k)-v(i,ay-2,k)))
          dv(i,ay-1,k) = sum(op%ar(1:6,ay-1)*(v(i,ay-5:ay,k)-v(i,ay-1,k)))+sum(op%ar(7:9,ay-1)*(vbr2y(i,1:3,k)-v(i,ay-1,k)))
          dv(i,ay  ,k) = sum(op%ar(1:5,ay  )*(v(i,ay-4:ay,k)-v(i,ay  ,k)))+sum(op%ar(6:9,ay  )*(vbr2y(i,1:4,k)-v(i,ay  ,k)))
        end do
      end do
      !$omp end target teams distribute parallel do
    endif
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do i=1,ax
        dv(i,1,k)=dv1(i,k)   ! supply lower solution
      end do; end do;
      !$omp end target teams distribute parallel do
    end if
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do i=1,ax
        dv(i,ay,k)=dv2(i,k)  ! supply upper solution
      end do; end do;
      !$omp end target teams distribute parallel do
    end if
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do
       !$omp end target teams distribute parallel do
      end if
      return
    endif
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do
       !$omp end target teams distribute parallel do
      end if
      return
    endif
    ! parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do k=1,az; do i=1,ax
        do j=1,2
          dvopy(j,i,k) = dv(i,j,k)
          dvopy(j+2,i,k) = dv(i,ay+j-2,k)
        end do
      end do; end do;
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do k=1,az; do i=1,ax
          dvopy(1:2,i,k) = zero
        end do; end do;
        !$omp end target teams distribute parallel do
      end if
      if( op%hi == MPI_PROC_NULL ) then
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do k=1,az; do i=1,ax
          dvopy(3:4,i,k) = zero
        end do; end do;
        !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopy)
      select case( dc )
      case( 1 ) ! mpi_allgather
        !$omp target data use_device_ptr(dvopy,dvoy) if(gpu_kernel==1)
        call mpi_allgather(dvopy,nsr,MPI_DOUBLE_PRECISION,dvoy,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvoy, np, ax, az )
        else
          call btrid_block4_lus( op%aa, dvoy, np, ax, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) then
          !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
          do k=1,az; do j=1,ay; do i=1,ax
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:2)*dvoy(3:4,i,k,op%lo)) 
          end do; end do; end do
          !$omp end target teams distribute parallel do
        end if
        if( op%hi /= MPI_PROC_NULL ) then 
          !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
          do k=1,az; do j=1,ay; do i=1,ax 
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,3:4)*dvoy(1:2,i,k,op%hi))
          end do; end do; end do
          !$omp end target teams distribute parallel do
        end if
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvopy,nsr,MPI_DOUBLE_PRECISION,dvoy,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvoy, np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvoy, np, ax, az )
          endif
        else
          dvoy = zero
        endif
        ! shuffle solution
        dvopy(3:4,:,:) = dvoy(1:2,:,:,0)
        do i=0,np-2
          dvopy(1:2,:,:) = dvoy(3:4,:,:,i)
          dvoy(3:4,:,:,i) = dvoy(1:2,:,:,i+1)
          dvoy(1:2,:,:,i+1) = dvopy(1:2,:,:)
        end do
        dvoy(1:2,:,:,0) = dvoy(3:4,:,:,np-1)
        dvoy(3:4,:,:,np-1) = dvopy(3:4,:,:)
        call mpi_scatter(dvoy,nsr,MPI_DOUBLE_PRECISION,dvopy,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvopy(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvopy(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvopy(:,i,k))
      end select
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do
        !$omp end target teams distribute parallel do
      end if
  end subroutine eval_compact_op1y_r4
  
  
  subroutine eval_compact_op1z_r4(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2) ! nor=4, nol=2, uses 1D op type
    implicit none
    integer, intent(in) :: ax,ay,az
    class(compact_op1_r4), intent(in) :: op
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double) :: vc
    integer :: nb,nsr,i,j,k,l
    integer :: nor,nir,nr,nol,nl,ni,np,dc  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) then
      	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
      	do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k) = zero
        end do; end do; end do;
        !$omp end target teams distribute parallel do
      end if
      if( op%null_option == 1 ) then
      	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
      	do k=1,az; do j=1,ay; do i=1,ax
      	  dv(i,j,k) = v(i,j,k)  ! filter / interp
      	end do; end do; end do;
        !$omp end target teams distribute parallel do
      end if
!      print *,'null op in z'
      RETURN
    END IF
    if( az /= op%m ) then
      print *,'*** error: mismatch in z operation size ***',az,op%m
      stop
    endif
    if( op%nor /= 4 ) then
      print *,'*** error: mismatch in z stencil size ***',4,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np  ; dc = op%directcom  
    ! explicit part
    ! ghost data
     if( np > 1 ) then  ! use parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay; do i=1,ax
        vbs2z(i,j,1:nor) = v(i,j,az-3:az)
        vbs1z(i,j,1:nor) = v(i,j,1:4)
      end do; end do;
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
      do j=1,ay; do i=1,ax
        vbr1z(i,j,1:nor) = v(i,j,az-3:az)
        vbr2z(i,j,1:nor) = v(i,j,1:4)
      end do; end do
      !$omp end target teams distribute parallel do
      if( debug ) print *,'periodic in z'
    endif

    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! antisymmetric BC, sumr != 0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,1) = sum(op%ar(5:9,1)*v(i,j,1:5))
          dv(i,j,2) = sum(op%ar(4:9,2)*v(i,j,1:6))
          dv(i,j,3) = sum(op%ar(3:9,3)*v(i,j,1:7))
          dv(i,j,4) = sum(op%ar(2:9,4)*v(i,j,1:8))
        end do
      end do
      !$omp end target teams distribute parallel do
     else ! symmetric/one-sided BC, assume sumr=0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,1) = sum(op%ar(6:9,1)*(v(i,j,2:5)-v(i,j,1)))
          dv(i,j,2) = sum(op%ar(5:9,2)*(v(i,j,2:6)-v(i,j,1)))
          dv(i,j,3) = sum(op%ar(4:9,3)*(v(i,j,2:7)-v(i,j,1)))
          dv(i,j,4) = sum(op%ar(3:9,4)*(v(i,j,2:8)-v(i,j,1)))
        end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,1) = sum(op%ar(1:4,1)*(vbr1z(i,j,1:4)-v(i,j,1)))+sum(op%ar(5:9,1)*(v(i,j,1:5)-v(i,j,1)))
          dv(i,j,2) = sum(op%ar(1:3,2)*(vbr1z(i,j,2:4)-v(i,j,2)))+sum(op%ar(4:9,2)*(v(i,j,1:6)-v(i,j,2)))
          dv(i,j,3) = sum(op%ar(1:2,3)*(vbr1z(i,j,3:4)-v(i,j,3)))+sum(op%ar(3:9,3)*(v(i,j,1:7)-v(i,j,3)))
          dv(i,j,4) = sum(op%ar(1:1,4)*(vbr1z(i,j,4:4)-v(i,j,4)))+sum(op%ar(2:9,4)*(v(i,j,1:8)-v(i,j,4)))
        end do
      end do
      !$omp end target teams distribute parallel do
    endif
    !$omp target teams distribute parallel do collapse(3) private(l,vc) if(gpu_kernel==1)
    do k=5,az-4
      do j=1,ay
        do i=1,ax
          dv(i,j,k) = zero
          vc = v(i,j,k)
          do l=1,9
            dv(i,j,k) = dv(i,j,k) + ( v(i,j,l+k-5) - vc )*op%ar(l,k)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! antisymmetric BC, sumr != 0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,az-3) = sum(op%ar(1:8,az-3)*v(i,j,az-7:az))
          dv(i,j,az-2) = sum(op%ar(1:7,az-2)*v(i,j,az-6:az))
          dv(i,j,az-1) = sum(op%ar(1:6,az-1)*v(i,j,az-5:az))
          dv(i,j,az  ) = sum(op%ar(1:5,az  )*v(i,j,az-4:az))
        end do
      end do
      !$omp end target teams distribute parallel do
     else ! symmetric/one-sided BC, assume sumr=0
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,az-3) = sum(op%ar(1:7,az-3)*(v(i,j,az-7:az-1)-v(i,j,az)))
          dv(i,j,az-2) = sum(op%ar(1:6,az-2)*(v(i,j,az-6:az-1)-v(i,j,az)))
          dv(i,j,az-1) = sum(op%ar(1:5,az-1)*(v(i,j,az-5:az-1)-v(i,j,az)))
          dv(i,j,az  ) = sum(op%ar(1:4,az  )*(v(i,j,az-4:az-1)-v(i,j,az)))
        end do
      end do
      !$omp end target teams distribute parallel do
     endif
    else
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay
        do i=1,ax
          dv(i,j,az-3) = sum(op%ar(1:8,az-3)*(v(i,j,az-7:az)-v(i,j,az-3)))+sum(op%ar(9:9,az-3)*(vbr2z(i,j,1:1)-v(i,j,az-3)))
          dv(i,j,az-2) = sum(op%ar(1:7,az-2)*(v(i,j,az-6:az)-v(i,j,az-2)))+sum(op%ar(8:9,az-2)*(vbr2z(i,j,1:2)-v(i,j,az-2)))
          dv(i,j,az-1) = sum(op%ar(1:6,az-1)*(v(i,j,az-5:az)-v(i,j,az-1)))+sum(op%ar(7:9,az-1)*(vbr2z(i,j,1:3)-v(i,j,az-1)))
          dv(i,j,az  ) = sum(op%ar(1:5,az  )*(v(i,j,az-4:az)-v(i,j,az  )))+sum(op%ar(6:9,az  )*(vbr2z(i,j,1:4)-v(i,j,az  )))
        end do
      end do
      !$omp end target teams distribute parallel do
    endif
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay; do i=1,ax
        dv(i,j,1)=dv1(i,j)   ! supply lower solution
      end do; end do
      !$omp end target teams distribute parallel do
    end if
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) then
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay; do i=1,ax
        dv(i,j,az)=dv2(i,j)  ! supply upper solution
      end do; end do
      !$omp end target teams distribute parallel do
    end if
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do;
        !$omp end target teams distribute parallel do
      end if
      return
    endif
   ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do;
        !$omp end target teams distribute parallel do
      end if
      return
    endif
    ! parallel solver
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay; do i=1,ax
        do k=1,2 
          dvopz(k,i,j) = dv(i,j,k)
          dvopz(k+2,i,j) = dv(i,j,az+k-2)
        end do
      end do; end do
      !$omp end target teams distribute parallel do
      if( op%lo == MPI_PROC_NULL ) then
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do j=1,ay; do i=1,ax
          dvopz(1:2,i,j) = zero
        end do; end do
        !$omp end target teams distribute parallel do
      end if
      if( op%hi == MPI_PROC_NULL ) then
        !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
        do j=1,ay; do i=1,ax
          dvopz(3:4,i,j) = zero
        end do; end do
        !$omp end target teams distribute parallel do
      end if
      nsr = size(dvopz)
      select case( dc )
      case( 1 ) ! mpi_allgather
        !$omp target data use_device_ptr(dvopz,dvoz) if(gpu_kernel==1)
        call mpi_allgather(dvopz,nsr,MPI_DOUBLE_PRECISION,dvoz,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        !$omp end target data
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvoz, np, ax, ay )
        else
          call btrid_block4_lus( op%aa, dvoz, np, ax, ay )
        endif
        if( op%lo /= MPI_PROC_NULL ) then
          !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
          do k=1,az; do j=1,ay; do i=1,ax
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:2)*dvoz(3:4,i,j,op%lo))
          end do; end do; end do;
          !$omp end target teams distribute parallel do
        end if
        if( op%hi /= MPI_PROC_NULL ) then 
          !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
          do k=1,az; do j=1,ay; do i=1,ax
            dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,3:4)*dvoz(1:2,i,j,op%hi))
          end do; end do; end do;
          !$omp end target teams distribute parallel do
        end if
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvopz,nsr,MPI_DOUBLE_PRECISION,dvoz,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvoz, np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvoz, np, ax, ay )
          endif
        else
          dvoz = zero
        endif
        ! shuffle solution
        dvopz(3:4,:,:) = dvoz(1:2,:,:,0)
        do i=0,np-2
          dvopz(1:2,:,:) = dvoz(3:4,:,:,i)
          dvoz(3:4,:,:,i) = dvoz(1:2,:,:,i+1)
          dvoz(1:2,:,:,i+1) = dvopz(1:2,:,:)
        end do
        dvoz(1:2,:,:,0) = dvoz(3:4,:,:,np-1)
        dvoz(3:4,:,:,np-1) = dvopz(3:4,:,:)
        call mpi_scatter(dvoz,nsr,MPI_DOUBLE_PRECISION,dvopz,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvopz(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvopz(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvopz(:,i,j))
      end select
      if( op%null_option == 1 ) then
        !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
        do k=1,az; do j=1,ay; do i=1,ax
          dv(i,j,k)=dv(i,j,k)+v(i,j,k)  ! filter / interp
        end do; end do; end do;
        !$omp end target teams distribute parallel do
      end if
  end subroutine eval_compact_op1z_r4


end module LES_compact_r4

