!===================================================================================================
! Copyright (c) 2018, Lawrence Livemore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! LLNL-CODE-749864
! This file is part of pyranda
! For details about use and distribution, please read: pyranda/LICENSE
!
! Written by: Britton J. Olson, olson45@llnl.gov
!===================================================================================================
module LES_compact
  USE iso_c_binding
  !USE MPI
  !USE LES_input, ONLY : bpp_lus_opt,use_ppent_opt,directcom,zerodx,zerody,zerodz
  USE LES_stencils
  USE LES_patch, ONLY : patch_type
  USE LES_comm,  ONLY : comm_type,mpierr,mpistatus,master
  use LES_pentadiagonal, ONLY : bpentLUD1, ppentLUD1, bpentLUS1, ppentlus, &
    bpentLUS3x, ppentLUS3x, bpentLUS3y, ppentLUS3y, bpentLUS3z, ppentLUS3z, &
    btrid_block4_lus, ptrid_block4_lus, btrid_block4_lud, ptrid_block4_lud
  use LES_pentadiagonal, only : bpentLUS2y
  USE LES_timers
  IMPLICIT NONE
  INCLUDE "mpif.h"
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double, one=1.0_c_double
  LOGICAL(c_bool) :: debug=.false.

  LOGICAL(c_bool), PARAMETER :: bpp_lus_opt=.true.   ! backward matrix compatibility
  LOGICAL(c_bool), PARAMETER :: use_ppent_opt=.true. ! backward matrix compatibility
  INTEGER(c_int), PARAMETER :: directcom = 1
  LOGICAL(c_bool), PARAMETER :: zerodx = .FALSE.
  LOGICAL(c_bool), PARAMETER :: zerody = .FALSE.
  LOGICAL(c_bool), PARAMETER :: zerodz = .FALSE.

  TYPE compact_op1  ! one direction, generic
    INTEGER(c_int) :: mo,no  ! operator size
    INTEGER(c_int) :: mv,nv  ! variable size
    INTEGER(c_int) :: m,n    ! mesh size
    INTEGER(c_int) :: nol,nor,nir,ncl,ncr,nci ! stencil
    INTEGER(c_int) :: nal,naa  ! parallel stencil
    INTEGER(c_int) :: directcom = 1  ! parallel solve
    INTEGER(c_int) :: null_option ! what to do if null_op (zero, copy, etc.)
    LOGICAL(c_bool) :: null_op    ! perform no op
    LOGICAL(c_bool) :: implicit_op  ! invert lhs
    LOGICAL(c_bool) :: periodic  ! periodicity
    INTEGER(c_int), dimension(2) :: bc,range  ! 1D bcs, range
    INTEGER(c_int) :: hash,np,id,lo,hi  ! 1D comm data
    real(KIND=c_double) :: d  ! nomimal grid spacing
    real(KIND=c_double) :: sh = zero ! shift from base grid
    real(kind=c_double), dimension(:,:), allocatable :: art ! transposed rhs (matrix compat.)
    real(kind=c_double), dimension(:,:), allocatable :: ar  ! rhs
    real(kind=c_double), dimension(:,:), allocatable :: al  ! lhs
    real(kind=c_double), dimension(:,:), allocatable :: rc  ! rhs overlap 
    real(kind=c_double), dimension(:,:,:,:), allocatable :: aa ! parallel weights
  contains   ! generic
    procedure :: setup => setup_compact_op1
    procedure :: remove => remove_compact_op1
    procedure :: evalx => eval_compact_op1x
    procedure :: evalbx => eval_compact_op1bx
    procedure :: ghostx => ghost_compact_op1x
    procedure :: evaly => eval_compact_op1y
    procedure :: evalz => eval_compact_op1z
  END TYPE compact_op1

  ! these extensions override the generic operations
  
  type, extends(compact_op1) :: compact_op1_d1  ! custom operator (backward compat.)
  contains
!    procedure :: evalx => eval_compact_op1x_d1t  ! matrix version (transposed rhs)
    procedure :: evalx => eval_compact_op1x_d1  ! matrix version
    procedure :: evaly => eval_compact_op1y_d1  ! matrix version
    procedure :: evalz => eval_compact_op1z_d1  ! matrix version
  end type
  
  type, extends(compact_op1) :: compact_op1_r3  ! custom operator (backward compat.)
  contains
    procedure :: evalx => eval_compact_op1x_r3  ! rhs stencil = 3
    procedure :: evaly => eval_compact_op1y_r3  ! rhs stencil = 3
    procedure :: evalz => eval_compact_op1z_r3  ! rhs stencil = 3
  end type
  
  type, extends(compact_op1) :: compact_op1_r4  ! custom operator (backward compat.)
  contains
    procedure :: evalx => eval_compact_op1x_r4  ! rhs stencil = 4
    procedure :: evaly => eval_compact_op1y_r4  ! rhs stencil = 4
    procedure :: evalz => eval_compact_op1z_r4  ! rhs stencil = 4
  end type
  
  type, extends(compact_op1) :: compact_cf1    ! custom operator
  contains
    procedure :: evalcfx => eval_compact_cf1x  ! different interface
    procedure :: evalcfy => eval_compact_cf1y  ! different interface
    procedure :: evalcfz => eval_compact_cf1z  ! different interface
  end type
  
  type, extends(compact_op1) :: compact_fc1    ! custom operator
  contains
    procedure :: evalfcx => eval_compact_fc1x  ! different interface
    procedure :: evalfcy => eval_compact_fc1y  ! different interface
    procedure :: evalfcz => eval_compact_fc1z  ! different interface
  end type
  
  TYPE control_type  ! from input or samrai
    LOGICAL(c_bool) :: null_opx=.false.,null_opy=.false.,null_opz=.false.    ! skip operations
    integer(c_int) :: d1spec=1,d2spec=1,d4spec=1,d8spec=1,gfspec=6,sfspec=1,tfspec=8  ! compact scheme
    integer(c_int) :: imcfspec=1,imfcspec=1,islspec=1,isrspec=2              ! compact scheme
    integer(c_int) :: amrcfspec=1,amrfcspec=2                                ! compact scheme
    INTEGER(c_int) :: directcom = 1  ! parallel solve
  END TYPE control_type
  
  type(control_type) :: control_data

  type compact_type  ! suite of operations
    type(control_type) :: control  ! controls
    integer(c_int) :: mbc(2,3,2), nop(3) ! (boundary, direction, symmetry), number of operators
    real(KIND=c_double) :: dx,dy,dz  ! nomimal grid spacing for derivatives
    type(compact_op1_d1), dimension(2) :: d1x,d1y,d1z
!    type(compact_op1_r3), dimension(2) :: d1x,d1y,d1z
    type(compact_op1_r3), dimension(2) :: d2x,d2y,d2z
    type(compact_op1_r3), dimension(2) :: d4x,d4y,d4z
    type(compact_op1_r4), dimension(2) :: d8x,d8y,d8z
    type(compact_op1_r4), dimension(2) :: gfx,gfy,gfz, sfx,sfy,sfz, tfx,tfy,tfz  ! gaussian, spectral, tophat filters
    type(compact_op1), dimension(2) :: isrx,isry,isrz,islx,isly,islz  ! arbitrarily right/left shifted interpolation
    type(compact_cf1), dimension(2) :: imcfx,imcfy,imcfz  ! center-face midpoint interpation N --> N+1
    type(compact_fc1), dimension(2) :: imfcx,imfcy,imfcz  ! face-center midpoint interpation N+1 --> N
    type(compact_op1), dimension(2) :: amrfcx,amrfcy,amrfcz, amrcfx,amrcfy,amrcfz ! invertible gaussian filters
  contains
    procedure :: setup => setup_compact_ops
    procedure :: remove => remove_compact_ops
  end type compact_type
  
  type(compact_type), pointer :: compact_ops => null()

!   TYPE line_type  ! 1D block and mesh data for compact setup
!     INTEGER(c_int)                 :: coordsys    ! coordinate system: 0=Cartesian, 1=Cylindrical, 2=Spherical, 3=Curvilinear
!     INTEGER(c_int)                 :: m,n         ! dimensions of mesh line
!     INTEGER(c_int)                 :: p           ! processor grid
!     REAL(c_double)                 :: d,bv1,bvn   ! bounding box (cell faces or nodes)
!     CHARACTER(KIND=c_char,LEN=4)   :: bc1,bcn     ! boundary conditions
!     INTEGER(c_int)                 :: isym        ! multipliers for symmetric and anti-symmetric boundaries
!     LOGICAL                        :: periodic    ! periodic directions
!     INTEGER         :: hash,np,id,lo,hi,range(2)  ! comm stuff -> [x,y,z]com,com_{np,id,lo,hi},range
!   END TYPE line_type
!   
!   type(line_type) :: xline,yline,zline

  TYPE mesh1_type  ! 1D
!    INTEGER(c_int) :: coordsys               ! coordinate system
    INTEGER(c_int) :: m,n                    ! mesh size  -> a[x,y,z], n[x,y,z]
    real(kind=c_double) :: d,bv1,bvn         ! grid spacing/metric -> d[x,y,z], [x,y,z][1,n]
    CHARACTER(KIND=c_char,LEN=4) :: bc1,bcn  ! boundary conditions (string) -> b[x,y,z][1,n]
  END TYPE mesh1_type
  
!   type(mesh1_type) :: xmesh_data,ymesh_data,zmesh_data

  TYPE comm1_type  ! 1D ! for testing
!    INTEGER(c_int) :: n,p  ! -> n[x,y,z], p[x,y,z]
    logical(c_bool) :: periodic  ! -> periodic[x,y,z]
    INTEGER(c_int) :: hash,np,id,lo,hi,range(2)  ! -> [x,y,z]com,com_{np,id,lo,hi},range
  END TYPE comm1_type

!   type(comm1_type) :: xcom_data,ycom_data,zcom_data

contains

  subroutine setup_compact_op1(op,wgt,com,msh,bc,null_op)  ! direction agnostic
    IMPLICIT NONE
    class(compact_op1), intent(out) :: op
    type(compact_weight), intent(in) :: wgt
    type(comm1_type), intent(in) :: com
    type(mesh1_type), intent(in) :: msh
    integer(c_int), intent(in) :: bc(2)
    logical(c_bool), intent(in) :: null_op
    real(kind=8), dimension(:,:), allocatable :: gar
    real(kind=8), dimension(:,:), allocatable :: gal
    real(kind=8), dimension(:,:), allocatable :: al1,al2
    real(kind=8), dimension(:,:), allocatable :: rop
    real(kind=8), dimension(:,:,:), allocatable :: ro
    integer :: m,n,np,nr,nl,ni,nol,nal,naa  ! local surrogates
    integer :: i,j
    ! copy weight data
    op%nol = wgt%nol ; op%ncl = wgt%ncl ; op%nci = wgt%nci
    op%nor = wgt%nor ; op%ncr = wgt%ncr ; op%nir = wgt%nir
    op%implicit_op = wgt%implicit_op ; op%sh = wgt%sh
    op%null_option = wgt%null_option
    nol = op%nol ; nl = op%ncl ; ni = op%nci ; nr = op%ncr
    ! copy 1D comm data
    op%periodic = com%periodic ; op%range = com%range
    op%hash = com%hash ; op%np = com%np ; op%id = com%id ; op%lo = com%lo ; op%hi = com%hi
    ! copy 1D mesh data
    op%m = msh%m ; op%n = msh%n ; op%mv = msh%m ; op%nv = msh%n ; op%d = msh%d
    select case( wgt%nst )  ! staggering
    case( 1 )    ! center to face
      if( op%hi == MPI_PROC_NULL ) op%m= op%m+1
      if( .not. op%periodic ) op%n = op%n+1
    case( -1 )   ! face to center
      op%mv = op%mv+1 ; op%nv = op%nv+1
    end select
    op%bc = bc
    if( op%n < 4 ) then
      op%null_op = .true.   ! skip by design
    else
      op%null_op = null_op  ! skip by request
    endif
    if( op%null_op ) return
    ! set up operator
    if( op%periodic ) then
      op%nal = nl+4
      op%naa = 4
    else
      op%nal = nl
      op%naa = 3
    endif
    nal = op%nal ; naa = op%naa ; m = op%m ; n = op%n ; np = op%np
    allocate( rop(ni,ni), al1(nol,nol), al2(nol,nol) )
    allocate( ro(ni,ni,0:np-1) )
    allocate( gal(nal,n), gar(nr,n) )
    allocate( op%ar(nr,m), op%al(m,nal) )
    allocate( op%art(m,nr) )
    gal = zero
    do i=1,n
      gal(1:nl,i) = wgt%ali
      gar(:,i) = wgt%ari
    end do
    if( .not. op%periodic ) then  ! patch USA options to ends
      gal(1:nl,1:wgt%nbc1) = wgt%alb1(:,:,bc(1))
      gal(1:nl,n-wgt%nbc2+1:n) = wgt%alb2(:,:,bc(2))
      gar(:,1:wgt%nbc1) = wgt%arb1(:,:,bc(1))
      gar(:,n-wgt%nbc2+1:n) = wgt%arb2(:,:,bc(2))
    endif
    ! *** debug
    if( debug ) then
    print *,'global coeff',m,n,nr,nal
    print *,wgt%ali
    print *,wgt%ari
    j = 1 ; if( op%periodic ) j = n
    do i=1,n,j
      print '(i6,19es14.6)',i,gal(:,i),gar(:,i)
    end do
    endif
    ! *** debug
    op%al = transpose(gal(:,op%range(1):op%range(1)+m-1))  ! lhs weights
    op%ar = gar(:,op%range(1):op%range(1)+m-1)             ! rhs weights
    op%art = transpose(op%ar)                   ! transposed rhs weights

    if( op%implicit_op ) then  ! lhs is needed
    al1 = zero
    al2 = zero
    if( op%lo /= MPI_PROC_NULL ) then ! interior corner 1
      do i=1,nol
        al1(i,i:nol) = op%al(i,1:nol-i+1)
      end do
      if( debug ) print *,'al1',al1
    endif
    if( op%hi /= MPI_PROC_NULL ) then ! interior corner 2
      do i=1,nol
        al2(nol-i+1,1:nol-i+1) = op%al(m-i+1,nl-nol+i:nl)
      end do
      if( debug ) print *,'al2',al2
    endif
    if( np == 1 ) then ! set up global LUD (m = n)
      if( op%periodic ) then
        call ppentLUD1(op%al,m)
      else
        call bpentLUD1(op%al,m)
      endif
    else  ! set up local LUD (m = n/np)
      allocate( op%rc(m,ni), op%aa(ni,ni,naa,0:np-1) )
      op%rc = zero
      call bpentLUD1(op%al,m)
      op%rc(1:nol,1:nol) = al1
      op%rc(m-nol+1:m,nol+1:ni) = al2
      do i=1,ni
        if( any(op%rc(:,i) /= zero) ) call bpentLUS1(op%al,op%rc(:,i),m)
        rop(1:nol,i) = op%rc(1:nol,i)
        rop(nol+1:ni,i) = op%rc(m-nol+1:m,i)
      end do
      if( op%lo == MPI_PROC_NULL ) rop(1:nol,:) = zero
      if( op%hi == MPI_PROC_NULL ) rop(nol+1:ni,:) = zero
      call startCOMM();call startCustom(1)
      call mpi_allgather(rop,ni*ni,MPI_DOUBLE_PRECISION,ro,ni*ni,MPI_DOUBLE_PRECISION,op%hash,mpierr)
      call endCOMM();call endCustom(1)
      op%aa = zero
      do j=0,np-1
        op%aa(:,nol+1:ni,1,j) = ro(:,1:nol,j)
        do i=1,ni
          op%aa(i,i,2,j) = one
        end do
        op%aa(:,1:nol,3,j) = ro(:,nol+1:ni,j)
      end do
      if( op%periodic ) then
        call ptrid_block4_lud( op%aa, np )
      else
        call btrid_block4_lud( op%aa, np )
      endif
    endif ! np
    endif ! implicit
    deallocate( gal, gar )
    deallocate( ro )
    deallocate( rop, al1, al2 )
  end subroutine setup_compact_op1

  subroutine remove_compact_op1(op)
    class(compact_op1), intent(out) :: op
    continue
  end subroutine remove_compact_op1

  function eval_compact_cf1x(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_cf1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1)+1,size(v,2),size(v,3)) :: dv
    real(kind=8), dimension(:,:), allocatable :: vbr,vbs
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    i = mx ; if( op%hi == MPI_PROC_NULL ) i=i+1
    if( i /= op%m ) then
      print *,'*** error: mismatch in operation size ***',i,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(nor,my,mz),vbr2(nor,my,mz) )
! ghost data
    allocate( vbs1(nor,my,mz),vbs2(nor,my,mz) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(mx-nor+1:mx,:,:)
      vbs1 = v(1:nor,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(mx-nor+1:mx,:,:)
      vbr2 = v(1:nor,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-nor+1:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:nor,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do j=1,my
      do i=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-i+1,i)*vbr1(i:nor,j,k))+sum(op%ar(nor-i+2:nr,i)*v(1:i+nir,j,k))
      end do
      do i=nor+1,mx-nor+1
        dv(i,j,k) = sum(op%ar(:,i)*v(i-nor:i+nir,j,k))
      end do
      do ii=1,nor
        i = mx-nor+ii
        dv(i+1,j,k) = sum(op%ar(1:nr-ii,i+1)*v(i-nir:mx,j,k))+sum(op%ar(nr-ii+1:nr,i+1)*vbr2(1:ii,j,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(mx+1,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
      dv(mx+1,:,:) = dv(1,:,:)
    else
      call bpentLUS3x(op%al,dv,size(op%al,1),size(dv,1),my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,my,mz), dvo(ni,my,mz,0:np-1) )
     forall(i=1:nol,j=1:my,k=1:mz) 
       dvop(i,j,k) = dv(i,j,k)
       dvop(i+nol,j,k) = dv(mx+i-nol,j,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, my, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, my, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:nol)*dvo(nol+1:ni,j,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,nol+1:ni)*dvo(1:nol,j,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, my, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, my, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
     end select
     deallocate( dvop, dvo )
     allocate( vbr(my,mz), vbs(my,mz) )
     vbs = dv(1,:,:)
     nsr = size(vbs)
     CALL startCOMM()
     call MPI_Sendrecv( vbs, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
     if( op%hi /= MPI_PROC_NULL ) dv(mx+1,:,:) = vbr
     deallocate( vbr, vbs )
  end function eval_compact_cf1x

  function eval_compact_cf1y(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_cf1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2)+1,size(v,3)) :: dv
    real(kind=8), dimension(:,:), allocatable :: vbr,vbs
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    j = my ; if( op%hi == MPI_PROC_NULL ) j=j+1
    if( j /= op%m ) then
      print *,'*** error: mismatch in operation size ***',j,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in y'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(mx,nor,mz),vbr2(mx,nor,mz) )
! ghost data
   allocate( vbs1(mx,nor,mz),vbs2(mx,nor,mz) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,mx-nor+1:mx,:)
     vbs1 = v(:,1:nor,:)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,my-nor+1:my,:)
      vbr2 = v(:,1:nor,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-nor+1:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:nor,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do i=1,mx
      do j=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-j+1,j)*vbr1(i,j:nor,k))+sum(op%ar(nor-j+2:nr,j)*v(i,1:j+nir,k))
      end do
      do j=nor+1,my-nor+1
        dv(i,j,k) = sum(op%ar(:,j)*v(i,j-nor:j+nir,k))
      end do
      do ii=1,nor
        j = my-nor+ii
        dv(i,j+1,k) = sum(op%ar(1:nr-ii,j+1)*v(i,j-nir:my,k))+sum(op%ar(nr-ii+1:nr,j+1)*vbr2(i,1:ii,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,my+1,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
      dv(:,my+1,:) = dv(:,1,:)
    else
      call bpentLUS3y(op%al,dv,size(op%al,1),mx,size(dv,2),mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,mz), dvo(ni,mx,mz,0:np-1) )
     forall(i=1:mx,j=1:nol,k=1:mz) 
       dvop(j,i,k) = dv(i,j,k)
       dvop(j+nol,i,k) = dv(i,my+j-nol,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:nol)*dvo(nol+1:ni,i,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,nol+1:ni)*dvo(1:nol,i,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
     end select
     deallocate( dvop, dvo )
     allocate( vbr(mx,mz), vbs(mx,mz) )
     vbs = dv(:,1,:)
     nsr = size(vbs)
     CALL startCOMM()
     call MPI_Sendrecv( vbs, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
     if( op%hi /= MPI_PROC_NULL ) dv(:,my+1,:) = vbr
     deallocate( vbr, vbs )
  end function eval_compact_cf1y

  function eval_compact_cf1z(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_cf1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)+1) :: dv
   real(kind=8), dimension(:,:), allocatable :: vbr,vbs
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    k = mz ; if( op%hi == MPI_PROC_NULL ) k=k+1
    if( k /= op%m ) then
      print *,'*** error: mismatch in operation size ***',k,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in z'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(mx,my,nor),vbr2(mx,my,nor) )
! ghost data
   allocate( vbs1(mx,my,nor),vbs2(mx,my,nor) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,:,mz-nor+1:mz)
     vbs1 = v(:,:,1:nor)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,mz-nor+1:mz)
      vbr2 = v(:,:,1:nor)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-nor+1:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:nor)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do j=1,my
    do i=1,mx
      do k=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-k+1,k)*vbr1(i,j,i:nor))+sum(op%ar(nor-k+2:nr,k)*v(i,j,1:k+nir))
      end do
      do k=nor+1,mz-nor+1
        dv(i,j,k) = sum(op%ar(:,k)*v(i,j,k-nor:k+nir))
      end do
      do ii=1,nor
        k = mz-nor+ii
        dv(i,j,k+1) = sum(op%ar(1:nr-ii,k+1)*v(i,j,k-nir:mz))+sum(op%ar(nr-ii+1:nr,k+1)*vbr2(i,j,1:ii))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,mz+1)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
      dv(:,:,mz+1) = dv(1,:,:)
   else
      call bpentLUS3z(op%al,dv,size(op%al,1),mx,my,size(dv,3))  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,my), dvo(ni,mx,my,0:np-1) )
     forall(i=1:mx,j=1:my,k=1:nol) 
       dvop(k,i,j) = dv(i,j,k)
       dvop(k+nol,i,j) = dv(i,j,mz+k-nol)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, my )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, my )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:nol)*dvo(nol+1:ni,i,j,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,nol+1:ni)*dvo(1:nol,i,j,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, my )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, my )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
     end select
     deallocate( dvop, dvo )
     allocate( vbr(mx,my), vbs(mx,my) )
     vbs = dv(:,:,1)
     nsr = size(vbs)
     CALL startCOMM()
     call MPI_Sendrecv( vbs, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
     if( op%hi /= MPI_PROC_NULL ) dv(:,:,mz+1) = vbr
     deallocate( vbr, vbs )
  end function eval_compact_cf1z

  function eval_compact_fc1x(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_fc1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1)-1,size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1)-1 ; my = size(v,2) ; mz = size(v,3)
    if( mx /= op%m ) then
      print *,'*** error: mismatch in operation size ***',mx,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(nor,my,mz),vbr2(nor,my,mz) )
! ghost data
   allocate( vbs1(nor,my,mz),vbs2(nor,my,mz) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(mx-nor+1:mx,:,:)
     vbs1 = v(2:nor+1,:,:)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(mx-nor+1:mx,:,:)
      vbr2 = v(2:nor+1,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-nor+1:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:nor,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do j=1,my
      do i=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-i+1,i)*vbr1(i:nor,j,k))+sum(op%ar(nor-i+2:nr,i)*v(1:i+nir,j,k))
      end do
      do i=nor+1,mx-nor
        dv(i,j,k) = sum(op%ar(:,i)*v(i-nor:i+nir,j,k))
      end do
      do ii=1,nor
        i = mx-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,i)*v(i-nir+1:mx+1,j,k))+sum(op%ar(nr-ii+1:nr,i)*vbr2(1:ii,j,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(mx-1,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
   else
      call bpentLUS3x(op%al,dv,size(op%al,1),mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,my,mz), dvo(ni,my,mz,0:np-1) )
     forall(i=1:nol,j=1:my,k=1:mz) 
       dvop(i,j,k) = dv(i,j,k)
       dvop(i+nol,j,k) = dv(mx+i-nol,j,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, my, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, my, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:nol)*dvo(nol+1:ni,j,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,nol+1:ni)*dvo(1:nol,j,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, my, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, my, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_fc1x

  function eval_compact_fc1y(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_fc1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2)-1,size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2)-1 ; mz = size(v,3)
    if( my /= op%m ) then
      print *,'*** error: mismatch in operation size ***',my,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in y'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(mx,nor,mz),vbr2(mx,nor,mz) )
! ghost data
   allocate( vbs1(mx,nor,mz),vbs2(mx,nor,mz) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,my-nor+1:my,:)
     vbs1 = v(:,2:nor+1,:)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,my-nor+1:my,:)
      vbr2 = v(:,2:nor+1,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-nor+1:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:nor,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do i=1,mx
      do j=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-j+1,j)*vbr1(i,j:nor,k))+sum(op%ar(nor-j+2:nr,j)*v(i,1:j+nir,k))
      end do
      do j=nor+1,my-nor
        dv(i,j,k) = sum(op%ar(:,j)*v(i,j-nor:j+nir,k))
      end do
      do ii=1,nor
        j = my-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,j)*v(i,j-nir+1:my+1,k))+sum(op%ar(nr-ii+1:nr,j)*vbr2(i,1:ii,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,my-1,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
   else
      call bpentLUS3y(op%al,dv,size(op%al,1),mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,mz), dvo(ni,mx,mz,0:np-1) )
     forall(i=1:mx,j=1:nol,k=1:mz) 
       dvop(j,i,k) = dv(i,j,k)
       dvop(j+nol,i,k) = dv(i,my+j-nol,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:nol)*dvo(nol+1:ni,i,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,nol+1:ni)*dvo(1:nol,i,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_fc1y

  function eval_compact_fc1z(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_fc1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)-1) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)-1
    if( mz /= op%m ) then
      print *,'*** error: mismatch in operation size ***',mz,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in z'
      RETURN
    END IF
    ! explicit part
    allocate( vbr1(mx,my,nor),vbr2(mx,my,nor) )
! ghost data
   allocate( vbs1(mx,my,nor),vbs2(mx,my,nor) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,:,mz-nor+1:mz)
     vbs1 = v(:,:,2:nor+1)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,mz-nor+1:mz)
      vbr2 = v(:,:,2:nor+1)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-nor+1:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:nor)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do j=1,my
    do i=1,mx
      do k=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-k+1,k)*vbr1(i,j,k:nor))+sum(op%ar(nor-k+2:nr,k)*v(i,j,1:k+nir))
      end do
      do k=nor+1,mz-nor
        dv(i,j,k) = sum(op%ar(:,k)*v(i,j,k-nor:k+nir))
      end do
      do ii=1,nor
        k = mz-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,k)*v(i,j,k-nir+1:mz+1))+sum(op%ar(nr-ii+1:nr,k)*vbr2(i,j,1:ii))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,mz-1)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
   else
      call bpentLUS3z(op%al,dv,size(op%al,1),mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,my), dvo(ni,mx,my,0:np-1) )
     forall(i=1:mx,j=1:my,k=1:nol) 
       dvop(k,i,j) = dv(i,j,k)
       dvop(k+nol,i,j) = dv(i,j,mz+k-nol)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, my )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, my )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:nol)*dvo(nol+1:ni,i,j,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,nol+1:ni)*dvo(1:nol,i,j,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, my )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, my )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_fc1z

  function eval_compact_op1x(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_op1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    if( mx /= op%m ) then
      print *,'*** error: mismatch in operation size ***',mx,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    ! explicit part
    allocate( vbr1(nor,my,mz),vbr2(nor,my,mz) )
! ghost data
   allocate( vbs1(nor,my,mz),vbs2(nor,my,mz) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(mx-nor+1:mx,:,:)
     vbs1 = v(1:nor,:,:)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(mx-nor+1:mx,:,:)
      vbr2 = v(1:nor,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-nor+1:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:nor,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do j=1,my
      do i=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-i+1,i)*vbr1(i:nor,j,k))+sum(op%ar(nor-i+2:nr,i)*v(1:i+nir,j,k))
      end do
      do i=nor+1,mx-nor
        dv(i,j,k) = sum(op%ar(:,i)*v(i-nor:i+nir,j,k))
      end do
      do ii=1,nor
        i = mx-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,i)*v(i-nir:mx,j,k))+sum(op%ar(nr-ii+1:nr,i)*vbr2(1:ii,j,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(mx,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,op%m,mx,my,mz)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,op%m,mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,my,mz), dvo(ni,my,mz,0:np-1) )
     forall(i=1:nol,j=1:my,k=1:mz) 
       dvop(i,j,k) = dv(i,j,k)
       dvop(i+nol,j,k) = dv(mx+i-nol,j,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, my, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, my, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:nol)*dvo(nol+1:ni,j,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,nol+1:ni)*dvo(1:nol,j,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, my, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, my, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_op1x

  function eval_compact_op1y(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_op1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in y'
      RETURN
    END IF
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    if( my /= op%m ) then
      print *,'*** error: mismatch in operation size ***',my,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    ! explicit part
    allocate( vbr1(mx,nor,mz),vbr2(mx,nor,mz) )
! ghost data
   allocate( vbs1(mx,nor,mz),vbs2(mx,nor,mz) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,my-nor+1:my,:)
     vbs1 = v(:,1:nor,:)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,my-nor+1:my,:)
      vbr2 = v(:,1:nor,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-nor+1:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:nor,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,mz
    do i=1,mx
      do j=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-j+1,j)*vbr1(i,j:nor,k))+sum(op%ar(nor-j+2:nr,j)*v(i,1:j+nir,k))
      end do
      do j=nor+1,my-nor
        dv(i,j,k) = sum(op%ar(:,j)*v(i,j-nor:j+nir,k))
      end do
      do ii=1,nor
        j = my-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,j)*v(i,j-nir:my,k))+sum(op%ar(nr-ii+1:nr,j)*vbr2(i,1:ii,k))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,my,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,op%m,mx,my,mz)  ! periodic solution on a single process
    else
      call bpentLUS3y(op%al,dv,op%m,mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,mz), dvo(ni,mx,mz,0:np-1) )
     forall(i=1:mx,j=1:nol,k=1:mz) 
       dvop(j,i,k) = dv(i,j,k)
       dvop(j+nol,i,k) = dv(i,my+j-nol,k)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:nol)*dvo(nol+1:ni,i,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,nol+1:ni)*dvo(1:nol,i,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_op1y

  function eval_compact_op1z(op,v,vb1,vb2,dv1,dv2) result(dv)  ! generalized, uses 1D op type
    implicit none
    class(compact_op1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in z'
      RETURN
    END IF
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    if( mz /= op%m ) then
      print *,'*** error: mismatch in operation size ***',mz,op%m
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np ;  id = op%id
    ! explicit part
    allocate( vbr1(mx,my,nor),vbr2(mx,my,nor) )
! ghost data
   allocate( vbs1(mx,my,nor),vbs2(mx,my,nor) )
   if( np > 1 ) then  ! use parallel solver
     vbs2 = v(:,:,mz-nor+1:mz)
     vbs1 = v(:,:,1:nor)
     nsr = size(vbs1)
     CALL startCOMM()
     call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                        vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                        op%hash, mpistatus, mpierr )
     call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                        vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                        op%hash, mpistatus, mpierr )
     CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,mz-nor+1:mz)
      vbr2 = v(:,:,1:nor)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-nor+1:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:nor)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do j=1,my
    do i=1,mx
      do k=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-k+1,k)*vbr1(i,j,k:nor))+sum(op%ar(nor-k+2:nr,k)*v(i,j,1:k+nir))
      end do
      do k=nor+1,mz-nor
        dv(i,j,k) = sum(op%ar(:,k)*v(i,j,k-nor:k+nir))
      end do
      do ii=1,nor
        k = mz-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,k)*v(i,j,k-nir:mz))+sum(op%ar(nr-ii+1:nr,k)*vbr2(i,j,1:ii))
      end do
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,mz)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,op%m,mx,my,mz)  ! periodic solution on a single process
    else
      call bpentLUS3z(op%al,dv,op%m,mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,mx,my), dvo(ni,mx,my,0:np-1) )
     forall(i=1:mx,j=1:my,k=1:nol) 
       dvop(k,i,j) = dv(i,j,k)
       dvop(k+nol,i,j) = dv(i,j,mz+k-nol)
     end forall
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, mx, my )
       else
         call btrid_block4_lus( op%aa, dvo, np, mx, my )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:nol)*dvo(nol+1:ni,i,j,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,nol+1:ni)*dvo(1:nol,i,j,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, mx, my )
         else
           call btrid_block4_lus( op%aa, dvo, np, mx, my )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_op1z

  function eval_compact_op1bx(op,v,vbr1,vbr2) result(dv)  ! generalized, uses 1D op type, ghost cells 
    implicit none
    class(compact_op1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in) :: vbr1,vbr2
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nsr,i,ii,j,k
    integer :: mx,my,mz,nor,nir,nr,nol,nl,ni,np,id  ! surrogates
    mx = size(v,1) ; my = size(v,2) ; mz = size(v,3)
    if( mx /= op%m .or. size(vbr1,1) /= op%nor .or. size(vbr2,1) /= op%nor ) then
      print *,'*** error: mismatch in operation size ***'
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np  ; id = op%id
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    ! explicit part
    do k=1,mz
    do j=1,my
      do i=1,nor
        dv(i,j,k) = sum(op%ar(1:nor-i+1,i)*vbr1(i:nor,j,k))+sum(op%ar(nor-i+2:nr,i)*v(1:i+nir,j,k))
      end do
      do i=nor+1,mx-nor
        dv(i,j,k) = sum(op%ar(:,i)*v(i-nor:i+nir,j,k))
      end do
      do ii=1,nor
        i = mx-nor+ii
        dv(i,j,k) = sum(op%ar(1:nr-ii,i)*v(i-nir:mx,j,k))+sum(op%ar(nr-ii+1:nr,i)*vbr2(1:ii,j,k))
      end do
    end do
    end do
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,size(op%al,1),mx,my,mz)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,size(op%al,1),mx,my,mz)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! use parallel solver
     allocate( dvop(ni,my,mz), dvo(ni,my,mz,0:np-1) )
     dvop(1:nol,:,:) = dv(1:nol,:,:)
     dvop(nol+1:ni,:,:) = dv(mx-nol+1:mx,:,:)
     if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
     if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
     nsr = size(dvop)
     select case( op%directcom )
     case( 1 ) ! mpi_allgather
        call startCOMM();call startCustom(1)
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        call endCOMM();call endCustom(1)
       if( op%periodic ) then
         call ptrid_block4_lus( op%aa, dvo, np, my, mz )
       else
         call btrid_block4_lus( op%aa, dvo, np, my, mz )
       endif
       if( op%lo /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:nol)*dvo(nol+1:ni,j,k,op%lo))
       if( op%hi /= MPI_PROC_NULL ) forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,nol+1:ni)*dvo(1:nol,j,k,op%hi))
    case( 2 ) ! mpi_gather/scatter
       call startCOMM();call startCustom(1)
       call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       call endCOMM();call endCustom(1)
       if( id == 0 ) then  ! only master solves
         if( op%periodic ) then
           call ptrid_block4_lus( op%aa, dvo, np, my, mz )
         else
           call btrid_block4_lus( op%aa, dvo, np, my, mz )
         endif
       else
         dvo = zero
       endif
       ! shuffle solution
       dvop(nol+1:ni,:,:) = dvo(1:nol,:,:,0)
       do i=0,np-2
         dvop(1:nol,:,:) = dvo(nol+1:ni,:,:,i)
         dvo(nol+1:ni,:,:,i) = dvo(1:nol,:,:,i+1)
         dvo(1:nol,:,:,i+1) = dvop(1:nol,:,:)
       end do
       dvo(1:nol,:,:,0) = dvo(nol+1:ni,:,:,np-1)
       dvo(nol+1:ni,:,:,np-1) = dvop(nol+1:ni,:,:)
       call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
       if( op%lo == MPI_PROC_NULL ) dvop(1:nol,:,:) = zero
       if( op%hi == MPI_PROC_NULL ) dvop(nol+1:ni,:,:) = zero
       forall(i=1:mx,j=1:my,k=1:mz) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
     end select
     deallocate( dvop, dvo )
  end function eval_compact_op1bx
  
  subroutine ghost_compact_op1x(op,v,vbr1,vbr2)  !,level,patch)
    implicit none
    class(compact_op1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), allocatable, intent(out) :: vbr1,vbr2
    real(kind=c_double), dimension(:,:,:), allocatable :: vbs1,vbs2
    integer :: nsr,nov
    nov = op%nor
    allocate( vbr1(nov,size(v,2),size(v,3)), vbr2(nov,size(v,2),size(v,3)) )
    if( op%np > 1 ) then  ! share data with neighbours
      allocate( vbs1(nov,size(v,2),size(v,3)), vbs2(nov,size(v,2),size(v,3)) )
      vbs2 = v(size(v,1)-nov+1:size(v,1),:,:)
      vbs1 = v(1:nov,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
      deallocate( vbs1,vbs2 )
    else if( op%periodic ) then
      vbr1 = v(size(v,1)-nov+1:size(v,1),:,:)
      vbr2 = v(1:nov,:,:)
    endif
    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data
  end subroutine ghost_compact_op1x

! "optimized" d1(t) operators for backward compatibility with matrix.f

  function eval_compact_op1x_d1t(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    real(kind=c_double), dimension(size(v,2),size(v,1)) :: dv_tran
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    character(len=160) :: filename
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
      dv = zero
!      print *,'null op in x'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
!---------------------------------------------------------------------------------------------------
    ! explicit part
    allocate( vbr1(3,ay,az),vbr2(3,ay,az) )
    allocate( vbs1(3,ay,az),vbs2(3,ay,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(ax-2:ax,:,:)
      vbs1 = v(1:3,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(ax-2:ax,:,:)
      vbr2 = v(1:3,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-2:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:3,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
!    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
!    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data

    if( op%lo == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do j=1,ay
      dv(1,j,k) = sum(op%art(1,4:7)*v(1:4,j,k))
      dv(2,j,k) = sum(op%art(2,3:7)*v(1:5,j,k))
      dv(3,j,k) = sum(op%art(3,2:7)*v(1:6,j,k))
     end do
     end do
     else
     do k=1,az
     do j=1,ay
      dv(1,j,k) = sum(op%art(1,5:7)*(v(2:4,j,k)-v(1,j,k)))
      dv(2,j,k) = sum(op%art(2,4:7)*(v(2:5,j,k)-v(1,j,k)))
      dv(3,j,k) = op%art(3,5)*(v(4,j,k)-v(2,j,k))+op%art(3,6)*(v(5,j,k)-v(1,j,k))+op%art(3,7)*(v(6,j,k)-v(1,j,k))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do j=1,ay
      dv(1,j,k) = op%art(1,5)*(v(2,j,k)-vbr1(3,j,k))+op%art(1,6)*(v(3,j,k)-vbr1(2,j,k))+op%art(1,7)*(v(4,j,k)-vbr1(1,j,k))
      dv(2,j,k) = op%art(2,5)*(v(3,j,k)-v(1,j,k))+op%art(2,6)*(v(4,j,k)-vbr1(3,j,k))+op%art(2,7)*(v(5,j,k)-vbr1(2,j,k))
      dv(3,j,k) = op%art(3,5)*(v(4,j,k)-v(2,j,k))+op%art(3,6)*(v(5,j,k)-v(1,j,k))+op%art(3,7)*(v(6,j,k)-vbr1(3,j,k))
    end do
     end do
    endif
    do k=1,az
    do j=1,ay
    do i=4,ax-3
      dv(i,j,k) = op%art(i,5)*(v(i+1,j,k)-v(i-1,j,k))+op%art(i,6)*(v(i+2,j,k)-v(i-2,j,k))+op%art(i,7)*(v(i+3,j,k)-v(i-3,j,k))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = sum(op%art(ax-2,1:6)*v(ax-5:ax,j,k))
      dv(ax-1,j,k) = sum(op%art(ax-1,1:5)*v(ax-4:ax,j,k))
      dv(ax,j,k)   = sum(op%art(ax,  1:4)*v(ax-3:ax,j,k))
     end do
     end do
     else
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = op%art(ax-2,1)*(v(ax-5,j,k)-v(ax,j,k))+op%art(ax-2,2)*(v(ax-4,j,k)-v(ax,j,k))+op%art(ax-2,3)*(v(ax-3,j,k)-v(ax-1,j,k))
      dv(ax-1,j,k) = sum(op%art(ax-1,1:4)*(v(ax-4:ax-1,j,k)-v(ax,j,k)))
      dv(ax,j,k)   = sum(op%art(ax,  1:3)*(v(ax-3:ax-1,j,k)-v(ax,j,k)))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = op%art(ax-2,5)*(v(ax-1,j,k)-v(ax-3,j,k))+op%art(ax-2,6)*(v(ax,j,k)-v(ax-4,j,k))+op%art(ax-2,7)*(vbr2(1,j,k)-v(ax-5,j,k))
      dv(ax-1,j,k) = op%art(ax-1,5)*(v(ax,j,k)-v(ax-2,j,k))+op%art(ax-1,6)*(vbr2(1,j,k)-v(ax-3,j,k))+op%art(ax-1,7)*(vbr2(2,j,k)-v(ax-4,j,k))
      dv(ax,j,k) = op%art(ax,5)*(vbr2(1,j,k)-v(ax-1,j,k))+op%art(ax,6)*(vbr2(2,j,k)-v(ax-2,j,k))+op%art(ax,7)*(vbr2(3,j,k)-v(ax-3,j,k))
     end do
     end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(ax,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      if (bpp_lus_opt) then
          do k = 1, az
             dv_tran = transpose(dv(:,:,k))
             call ppentlus_f77(ay,ax,1,op%al(1,1),dv_tran(1,1))
             dv(:,:,k) = transpose(dv_tran)
          end do
      else
         call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
      end if
    else
      if (bpp_lus_opt) then
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
      endif
    endif
    if( np > 1 ) then  ! use parallel solver
      allocate( dvop(4,ay,az), dvo(4,ay,az,0:np-1) )
      dvop(1:2,:,:) = dv(1:2,:,:)
      dvop(3:4,:,:) = dv(ax-1:ax,:,:)
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ay, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ay, az )
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvo(3,j,k,op%lo) - op%rc(i,2)*dvo(4,j,k,op%lo) &
                        &           - op%rc(i,3)*dvo(1,j,k,op%hi) - op%rc(i,4)*dvo(2,j,k,op%hi)
           end do
           end do
           end do
        else if (op%lo /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvo(3,j,k,op%lo) - op%rc(i,2)*dvo(4,j,k,op%lo)
           end do
           end do
           end do
        else if (op%hi /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,3)*dvo(1,j,k,op%hi) - op%rc(i,4)*dvo(2,j,k,op%hi)
           end do
           end do
           end do
        end if
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ay, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
      end select
      deallocate( dvop, dvo )
    endif
  end function eval_compact_op1x_d1t

  function eval_compact_op1x_d1(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    real(kind=c_double), dimension(size(v,2),size(v,1)) :: dv_tran
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    character(len=160) :: filename
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
      dv = zero
!      print *,'null op in x'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
!---------------------------------------------------------------------------------------------------
    ! explicit part
    allocate( vbr1(3,ay,az),vbr2(3,ay,az) )
    allocate( vbs1(3,ay,az),vbs2(3,ay,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(ax-2:ax,:,:)
      vbs1 = v(1:3,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(ax-2:ax,:,:)
      vbr2 = v(1:3,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-2:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:3,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
!    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
!    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data
     if( op%lo == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do j=1,ay
      dv(1,j,k) = sum(op%ar(4:7,1)*v(1:4,j,k))
      dv(2,j,k) = sum(op%ar(3:7,2)*v(1:5,j,k))
      dv(3,j,k) = sum(op%ar(2:7,3)*v(1:6,j,k))
     end do
     end do
     else
     do k=1,az
     do j=1,ay
      dv(1,j,k) = sum(op%ar(5:7,1)*(v(2:4,j,k)-v(1,j,k)))
      dv(2,j,k) = sum(op%ar(4:7,2)*(v(2:5,j,k)-v(1,j,k)))
      dv(3,j,k) = op%ar(5,3)*(v(4,j,k)-v(2,j,k))+op%ar(6,3)*(v(5,j,k)-v(1,j,k))+op%ar(7,3)*(v(6,j,k)-v(1,j,k))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do j=1,ay
      dv(1,j,k) = op%ar(5,1)*(v(2,j,k)-vbr1(3,j,k))+op%ar(6,1)*(v(3,j,k)-vbr1(2,j,k))+op%ar(7,1)*(v(4,j,k)-vbr1(1,j,k))
      dv(2,j,k) = op%ar(5,2)*(v(3,j,k)-v(1,j,k))+op%ar(6,2)*(v(4,j,k)-vbr1(3,j,k))+op%ar(7,2)*(v(5,j,k)-vbr1(2,j,k))
      dv(3,j,k) = op%ar(5,3)*(v(4,j,k)-v(2,j,k))+op%ar(6,3)*(v(5,j,k)-v(1,j,k))+op%ar(7,3)*(v(6,j,k)-vbr1(3,j,k))
     end do
     end do
    endif
    do k=1,az
    do j=1,ay
    do i=4,ax-3
      dv(i,j,k) = op%ar(5,i)*(v(i+1,j,k)-v(i-1,j,k))+op%ar(6,i)*(v(i+2,j,k)-v(i-2,j,k))+op%ar(7,i)*(v(i+3,j,k)-v(i-3,j,k))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL ) then ! boundary weights
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*v(ax-5:ax,j,k))
      dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*v(ax-4:ax,j,k))
      dv(ax,j,k)   = sum(op%ar(1:4,ax  )*v(ax-3:ax,j,k))
     end do
     end do
     else
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = op%ar(1,ax-2)*(v(ax-5,j,k)-v(ax,j,k))+op%ar(2,ax-2)*(v(ax-4,j,k)-v(ax,j,k))+op%ar(3,ax-2)*(v(ax-3,j,k)-v(ax-1,j,k))
      dv(ax-1,j,k) = sum(op%ar(1:4,ax-1)*(v(ax-4:ax-1,j,k)-v(ax,j,k)))
      dv(ax,j,k)   = sum(op%ar(1:3,ax  )*(v(ax-3:ax-1,j,k)-v(ax,j,k)))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do j=1,ay
      dv(ax-2,j,k) = op%ar(5,ax-2)*(v(ax-1,j,k)-v(ax-3,j,k))+op%ar(6,ax-2)*(v(ax,j,k)-v(ax-4,j,k))+op%ar(7,ax-2)*(vbr2(1,j,k)-v(ax-5,j,k))
      dv(ax-1,j,k) = op%ar(5,ax-1)*(v(ax,j,k)-v(ax-2,j,k))+op%ar(6,ax-1)*(vbr2(1,j,k)-v(ax-3,j,k))+op%ar(7,ax-1)*(vbr2(2,j,k)-v(ax-4,j,k))
      dv(ax,j,k)   = op%ar(5,ax  )*(vbr2(1,j,k)-v(ax-1,j,k))+op%ar(6,ax)*(vbr2(2,j,k)-v(ax-2,j,k))+op%ar(7,ax  )*(vbr2(3,j,k)-v(ax-3,j,k))
     end do
     end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(ax,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      if (bpp_lus_opt) then
          do k = 1, az
             dv_tran = transpose(dv(:,:,k))
             call ppentlus_f77(ay,ax,1,op%al(1,1),dv_tran(1,1))
             dv(:,:,k) = transpose(dv_tran)
          end do
      else
         call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
      end if
    else
      if (bpp_lus_opt) then
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
      endif
    endif
    if( np > 1 ) then  ! use parallel solver
      allocate( dvop(4,ay,az), dvo(4,ay,az,0:np-1) )
      dvop(1:2,:,:) = dv(1:2,:,:)
      dvop(3:4,:,:) = dv(ax-1:ax,:,:)
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ay, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ay, az )
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvo(3,j,k,op%lo) - op%rc(i,2)*dvo(4,j,k,op%lo) &
                        &           - op%rc(i,3)*dvo(1,j,k,op%hi) - op%rc(i,4)*dvo(2,j,k,op%hi)
           end do
           end do
           end do
        else if (op%lo /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,1)*dvo(3,j,k,op%lo) - op%rc(i,2)*dvo(4,j,k,op%lo)
           end do
           end do
           end do
        else if (op%hi /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(i,3)*dvo(1,j,k,op%hi) - op%rc(i,4)*dvo(2,j,k,op%hi)
           end do
           end do
           end do
        end if
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ay, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
      end select
      deallocate( dvop, dvo )
    endif
  end function eval_compact_op1x_d1

  function eval_compact_op1y_d1(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
      dv = zero
!      print *,'null op in y'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    if( ay /= op%m ) then
      print *,'*** error: mismatch in y operation size ***',ax,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in y stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
!---------------------------------------------------------------------------------------------------      
! ghost data
    allocate( vbr1(ax,3,az),vbr2(ax,3,az) )
    allocate( vbs1(ax,3,az),vbs2(ax,3,az) )
    ! explicit part
    if( op%np > 1 ) then  ! use parallel solver
      vbs2 = v(:,ay-2:ay,:)
      vbs1 = v(:,1:3,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,ay-2:ay,:)
      vbr2 = v(:,1:3,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-2:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:3,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
!    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
!    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data
    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do i=1,ax
      dv(i,1,k) = sum(op%ar(4:7,1)*v(i,1:4,k))
      dv(i,2,k) = sum(op%ar(3:7,2)*v(i,1:5,k))
      dv(i,3,k) = sum(op%ar(2:7,3)*v(i,1:6,k))
     end do
     end do
     else
     do k=1,az
     do i=1,ax
      dv(i,1,k) = sum(op%ar(5:7,1)*(v(i,2:4,k)-v(i,1,k)))
      dv(i,2,k) = sum(op%ar(4:7,2)*(v(i,2:5,k)-v(i,1,k)))
      dv(i,3,k) = op%ar(5,3)*(v(i,4,k)-v(i,2,k))+op%ar(6,3)*(v(i,5,k)-v(i,1,k))+op%ar(7,3)*(v(i,6,k)-v(i,1,k))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do i=1,ax
      dv(i,1,k) = op%ar(5,1)*(v(i,2,k)-vbr1(i,3,k))+op%ar(6,1)*(v(i,3,k)-vbr1(i,2,k))+op%ar(7,1)*(v(i,4,k)-vbr1(i,1,k))
      dv(i,2,k) = op%ar(5,2)*(v(i,3,k)-v(i,1,k))+op%ar(6,2)*(v(i,4,k)-vbr1(i,3,k))+op%ar(7,2)*(v(i,5,k)-vbr1(i,2,k))
      dv(i,3,k) = op%ar(5,3)*(v(i,4,k)-v(i,2,k))+op%ar(6,3)*(v(i,5,k)-v(i,1,k))+op%ar(7,3)*(v(i,6,k)-vbr1(i,3,k))
     end do
     end do
    endif
    do k=1,az
    do j=4,ay-3
    do i=1,ax
      dv(i,j,k) = op%ar(5,j)*(v(i,j+1,k)-v(i,j-1,k))+op%ar(6,j)*(v(i,j+2,k)-v(i,j-2,k))+op%ar(7,j)*(v(i,j+3,k)-v(i,j-3,k))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
     do k=1,az
     do i=1,ax
      dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*v(i,ay-5:ay,k))
      dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*v(i,ay-4:ay,k))
      dv(i,ay,k)   = sum(op%ar(1:4,ay  )*v(i,ay-3:ay,k))
     end do
     end do
     else
     do k=1,az
     do i=1,ax
      dv(i,ay-2,k) = op%ar(1,ay-2)*(v(i,ay-5,k)-v(i,ay,k))+op%ar(2,ay-2)*(v(i,ay-4,k)-v(i,ay,k))+op%ar(3,ay-2)*(v(i,ay-3,k)-v(i,ay-1,k))
      dv(i,ay-1,k) = sum(op%ar(1:4,ay-1)*(v(i,ay-4:ay-1,k)-v(i,ay,k)))
      dv(i,ay,k)   = sum(op%ar(1:3,ay  )*(v(i,ay-3:ay-1,k)-v(i,ay,k)))
     end do
     end do
     endif
    else ! centered interior weights
     do k=1,az
     do i=1,ax
      dv(i,ay-2,k) = op%ar(5,ay-2)*(v(i,ay-1,k)-v(i,ay-3,k))+op%ar(6,ay-2)*(v(i,ay,k)-v(i,ay-4,k))+op%ar(7,ay-2)*(vbr2(i,1,k)-v(i,ay-5,k))
      dv(i,ay-1,k) = op%ar(5,ay-1)*(v(i,ay,k)-v(i,ay-2,k))+op%ar(6,ay-1)*(vbr2(i,1,k)-v(i,ay-3,k))+op%ar(7,ay-1)*(vbr2(i,2,k)-v(i,ay-4,k))
      dv(i,ay,k)   = op%ar(5,ay  )*(vbr2(i,1,k)-v(i,ay-1,k))+op%ar(6,ay)*(vbr2(i,2,k)-v(i,ay-2,k))+op%ar(7,ay  )*(vbr2(i,3,k)-v(i,ay-3,k))
     end do
     end do
    endif
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1     ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,ay+1,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( op%np == 1 .and. op%periodic ) then
       if (bpp_lus_opt) then
          call ppentlus_f77(ax,ay,az,op%al(1,1),dv(1,1,1))  ! periodic solution on a single process
       else
          call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
       end if
    else
       call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( op%np > 1 ) then  ! use parallel solver
      allocate( dvop(4,ax,az), dvo(4,ax,az,0:op%np-1) )
      do k=1,az
      do i=1,ax
        dvop(1:2,i,k) = dv(i,1:2,k)
        dvop(3:4,i,k) = dv(i,ay-1:ay,k)
      end do
      end do
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, op%np, ax, az )
        else
          call btrid_block4_lus( op%aa, dvo, op%np, ax, az )
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(j,1)*dvo(3,i,k,op%lo) - op%rc(j,2)*dvo(4,i,k,op%lo) &
                        &           - op%rc(j,3)*dvo(1,i,k,op%hi) - op%rc(j,4)*dvo(2,i,k,op%hi)
           end do
           end do
           end do
        else if (op%lo /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(j,1)*dvo(3,i,k,op%lo) - op%rc(j,2)*dvo(4,i,k,op%lo)
           end do
           end do
           end do
        else if (op%hi /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(j,3)*dvo(1,i,k,op%hi) - op%rc(j,4)*dvo(2,i,k,op%hi)
           end do
           end do
           end do
        end if
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, op%np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvo, op%np, ax, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,op%np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,op%np-1)
        dvo(3:4,:,:,op%np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
      end select
      deallocate( dvop, dvo )
    endif
  end function eval_compact_op1y_d1

  function eval_compact_op1z_d1(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_d1), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
!---------------------------------------------------------------------------------------------------
    IF (op%null_op) THEN
      dv = zero
!      print *,'null op in z'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
!---------------------------------------------------------------------------------------------------      
    ! explicit part
! ghost data
    allocate( vbr1(ax,ay,3),vbr2(ax,ay,3) )
    allocate( vbs1(ax,ay,3),vbs2(ax,ay,3) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,:,az-2:az)
      vbs1 = v(:,:,1:3)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,az-2:az)
      vbr2 = v(:,:,1:3)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-2:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:3)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
!    if( op%lo == MPI_PROC_NULL ) vbr1 = zero  ! lower non-periodic boundary, no data
!    if( op%hi == MPI_PROC_NULL ) vbr2 = zero  ! upper non-periodic boundary, no data
    if( op%lo == MPI_PROC_NULL ) then
     if( op%bc(1) == -1 ) then ! this BC is antisymmetric
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = sum(op%ar(4:7,1)*v(i,j,1:4))
      dv(i,j,2) = sum(op%ar(3:7,2)*v(i,j,1:5))
      dv(i,j,3) = sum(op%ar(2:7,3)*v(i,j,1:6))
     end do
     end do
     else
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = sum(op%ar(5:7,1)*(v(i,j,2:4)-v(i,j,1)))
      dv(i,j,2) = sum(op%ar(4:7,2)*(v(i,j,2:5)-v(i,j,1)))
      dv(i,j,3) = op%ar(5,3)*(v(i,j,4)-v(i,j,2))+op%ar(6,3)*(v(i,j,5)-v(i,j,1))+op%ar(7,3)*(v(i,j,6)-v(i,j,1))
     end do
     end do
     endif
    else ! centered interior weights
     do j=1,ay
     do i=1,ax
      dv(i,j,1) = op%ar(5,1)*(v(i,j,2)-vbr1(i,j,3))+op%ar(6,1)*(v(i,j,3)-vbr1(i,j,2))+op%ar(7,1)*(v(i,j,4)-vbr1(i,j,1))
      dv(i,j,2) = op%ar(5,2)*(v(i,j,3)-v(i,j,1))+op%ar(6,2)*(v(i,j,4)-vbr1(i,j,3))+op%ar(7,2)*(v(i,j,5)-vbr1(i,j,2))
      dv(i,j,3) = op%ar(5,3)*(v(i,j,4)-v(i,j,2))+op%ar(6,3)*(v(i,j,5)-v(i,j,1))+op%ar(7,3)*(v(i,j,6)-vbr1(i,j,3))
     end do
     end do
    endif
    do k=4,az-3
    do j=1,ay
    do i=1,ax
      dv(i,j,k) = op%ar(5,k)*(v(i,j,k+1)-v(i,j,k-1))+op%ar(6,k)*(v(i,j,k+2)-v(i,j,k-2))+op%ar(7,k)*(v(i,j,k+3)-v(i,j,k-3))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL ) then
     if( op%bc(2) == -1 ) then ! this BC is antisymmetric
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = sum(op%ar(1:6,az-2)*v(i,j,az-5:az))
      dv(i,j,az-1) = sum(op%ar(1:5,az-1)*v(i,j,az-4:az))
      dv(i,j,az)   = sum(op%ar(1:4,az  )*v(i,j,az-3:az))
     end do
     end do
     else
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = op%ar(1,az-2)*(v(i,j,az-5)-v(i,j,az))+op%ar(2,az-2)*(v(i,j,az-4)-v(i,j,az))+op%ar(3,az-2)*(v(i,j,az-3)-v(i,j,az-1))
      dv(i,j,az-1) = sum(op%ar(1:4,az-1)*(v(i,j,az-4:az-1)-v(i,j,az)))
      dv(i,j,az)   = sum(op%ar(1:3,az  )*(v(i,j,az-3:az-1)-v(i,j,az)))
     end do
     end do
     endif
    else ! centered interior weights
     do j=1,ay
     do i=1,ax
      dv(i,j,az-2) = op%ar(5,az-2)*(v(i,j,az-1)-v(i,j,az-3))+op%ar(6,az-2)*(v(i,j,az)-v(i,j,az-4))+op%ar(7,az-2)*(vbr2(i,j,1)-v(i,j,az-5))
      dv(i,j,az-1) = op%ar(5,az-1)*(v(i,j,az)-v(i,j,az-2))+op%ar(6,az-1)*(vbr2(i,j,1)-v(i,j,az-3))+op%ar(7,az-1)*(vbr2(i,j,2)-v(i,j,az-4))
      dv(i,j,az)   = op%ar(5,az  )*(vbr2(i,j,1)-v(i,j,az-1))+op%ar(6,az)*(vbr2(i,j,2)-v(i,j,az-2))+op%ar(7,az  )*(vbr2(i,j,3)-v(i,j,az-3))
     end do
     end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,az)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( op%np == 1 .and. op%periodic ) then
       if (bpp_lus_opt) then
          call ppentlus(3,use_ppent_opt,op%al,dv)
       else
          call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
       end if
    else
      call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( op%np > 1 ) then  ! use parallel solver
      allocate( dvop(4,ax,ay), dvo(4,ax,ay,0:op%np-1) )
      do k=1,2
        dvop(k,:,:) = dv(:,:,k)
        dvop(2+k,:,:) = dv(:,:,az-2+k)
      end do
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, op%np, ax, ay )
        else
          call btrid_block4_lus( op%aa, dvo, op%np, ax, ay )
        endif
        if ((op%lo /= MPI_PROC_NULL) .and. (op%hi /= MPI_PROC_NULL)) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,1)*dvo(3,i,j,op%lo) - op%rc(k,2)*dvo(4,i,j,op%lo) &
                        &           - op%rc(k,3)*dvo(1,i,j,op%hi) - op%rc(k,4)*dvo(2,i,j,op%hi)
           end do
           end do
           end do
        else if (op%lo /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,1)*dvo(3,i,j,op%lo) - op%rc(k,2)*dvo(4,i,j,op%lo)
           end do
           end do
           end do
        else if (op%hi /= MPI_PROC_NULL) then
           do k = 1, az
           do j = 1, ay
           do i = 1, ax
              dv(i,j,k) = dv(i,j,k) - op%rc(k,3)*dvo(1,i,j,op%hi) - op%rc(k,4)*dvo(2,i,j,op%hi)
           end do
           end do
           end do
        end if
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, op%np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvo, op%np, ax, ay )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,op%np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,op%np-1)
        dvo(3:4,:,:,op%np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
      end select
      deallocate( dvop, dvo )
    endif
  end function eval_compact_op1z_d1

! "optimized" r3 operators for backward compatibility with matrix.f

  function eval_compact_op1x_r3(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
    ! explicit part
    allocate( vbr1(3,ay,az),vbr2(3,ay,az) )
! ghost data
   allocate( vbs1(3,ay,az),vbs2(3,ay,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(ax-2:ax,:,:)
      vbs1 = v(1:3,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(ax-2:ax,:,:)
      vbr2 = v(1:3,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-2:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:3,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,az
    do j=1,ay
      dv(1,j,k) = sum(op%ar(1:3,1)*vbr1(1:3,j,k))+sum(op%ar(4:7,1)*v(1:4,j,k))
      dv(2,j,k) = sum(op%ar(1:2,2)*vbr1(2:3,j,k))+sum(op%ar(3:7,2)*v(1:5,j,k))
      dv(3,j,k) = sum(op%ar(1:1,3)*vbr1(3:3,j,k))+sum(op%ar(2:7,3)*v(1:6,j,k))
    end do
    end do
    do k=1,az
    do j=1,ay
      do i=4,ax-3
        dv(i,j,k) = sum(op%ar(:,i)*v(i-3:i+3,j,k))
      end do
    end do
    end do
    do k=1,az
    do j=1,ay
      dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*v(ax-5:ax,j,k))+sum(op%ar(7:7,ax-2)*vbr2(1:1,j,k))
      dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*v(ax-4:ax,j,k))+sum(op%ar(6:7,ax-1)*vbr2(1:2,j,k))
      dv(ax  ,j,k) = sum(op%ar(1:4,ax  )*v(ax-3:ax,j,k))+sum(op%ar(5:7,ax  )*vbr2(1:3,j,k))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(ax,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ay,az), dvo(4,ay,az,0:np-1) )
      forall(i=1:2,j=1:ay,k=1:az) 
        dvop(i,j,k) = dv(i,j,k)
        dvop(i+2,j,k) = dv(ax+i-2,j,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ay, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ay, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:2)*dvo(3:4,j,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,3:4)*dvo(1:2,j,k,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ay, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1x_r3

  function eval_compact_op1y_r3(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in y'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    if( ay /= op%m ) then
      print *,'*** error: mismatch in y operation size ***',ay,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in y stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
    ! explicit part
    allocate( vbr1(ax,3,az),vbr2(ax,3,az) )
! ghost data
   allocate( vbs1(ax,3,az),vbs2(ax,3,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,ay-2:ay,:)
      vbs1 = v(:,1:3,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,ay-2:ay,:)
      vbr2 = v(:,1:3,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-2:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:3,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,az
    do i=1,ax
      dv(i,1,k) = sum(op%ar(1:3,1)*vbr1(i,1:3,k))+sum(op%ar(4:7,1)*v(i,1:4,k))
      dv(i,2,k) = sum(op%ar(1:2,2)*vbr1(i,2:3,k))+sum(op%ar(3:7,2)*v(i,1:5,k))
      dv(i,3,k) = sum(op%ar(1:1,3)*vbr1(i,3:3,k))+sum(op%ar(2:7,3)*v(i,1:6,k))
    end do
    end do
    do k=1,az
    do j=4,ay-3
    do i=1,ax
        dv(i,j,k) = sum(op%ar(:,j)*v(i,j-3:j+3,k))
    end do
    end do
    end do
    do k=1,az
    do i=1,ax
      dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*v(i,ay-5:ay,k))+sum(op%ar(7:7,ay-2)*vbr2(i,1:1,k))
      dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*v(i,ay-4:ay,k))+sum(op%ar(6:7,ay-1)*vbr2(i,1:2,k))
      dv(i,ay  ,k) = sum(op%ar(1:4,ay  )*v(i,ay-3:ay,k))+sum(op%ar(5:7,ay  )*vbr2(i,1:3,k))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,ay,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ax,az), dvo(4,ax,az,0:np-1) )
      forall(i=1:ax,j=1:2,k=1:az) 
        dvop(j,i,k) = dv(i,j,k)
        dvop(j+2,i,k) = dv(i,ay+j-2,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:2)*dvo(3:4,i,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,3:4)*dvo(1:2,i,k,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1y_r3

  function eval_compact_op1z_r3(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in z'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    if( az /= op%m ) then
      print *,'*** error: mismatch in z operation size ***',az,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in z stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
    ! explicit part
! ghost data
    allocate( vbr1(ax,ay,3),vbr2(ax,ay,3) )
    allocate( vbs1(ax,ay,3),vbs2(ax,ay,3) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,:,az-2:az)
      vbs1 = v(:,:,1:3)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,az-2:az)
      vbr2 = v(:,:,1:3)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-2:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:3)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do j=1,ay
    do i=1,ax
      dv(i,j,1) = sum(op%ar(1:3,1)*vbr1(i,j,1:3))+sum(op%ar(4:7,1)*v(i,j,1:4))
      dv(i,j,2) = sum(op%ar(1:2,2)*vbr1(i,j,2:3))+sum(op%ar(3:7,2)*v(i,j,1:5))
      dv(i,j,3) = sum(op%ar(1:1,3)*vbr1(i,j,3:3))+sum(op%ar(2:7,3)*v(i,j,1:6))
    end do
    end do
    do k=4,az-3
    do j=1,ay
    do i=1,ax
        dv(i,j,k) = sum(op%ar(:,k)*v(i,j,k-3:k+3))
    end do
    end do
    end do
    do j=1,ay
    do i=1,ax
      dv(i,j,az-2) = sum(op%ar(1:6,az-2)*v(i,j,az-5:az))+sum(op%ar(7:7,az-2)*vbr2(i,j,1:1))
      dv(i,j,az-1) = sum(op%ar(1:5,az-1)*v(i,j,az-4:az))+sum(op%ar(6:7,az-1)*vbr2(i,j,1:2))
      dv(i,j,az  ) = sum(op%ar(1:4,az  )*v(i,j,az-3:az))+sum(op%ar(5:7,az  )*vbr2(i,j,1:3))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,az)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ax,ay), dvo(4,ax,ay,0:np-1) )
      forall(i=1:ax,j=1:ay,k=1:2) 
        dvop(k,i,j) = dv(i,j,k)
        dvop(k+2,i,j) = dv(i,j,az+k-2)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, ay )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:2)*dvo(3:4,i,j,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,3:4)*dvo(1:2,i,j,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, ay )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1z_r3

! "optimized" r4 operators for backward compatibility with matrix.f

  function eval_compact_op1x_r4(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=4, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r4), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in x'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
    ! explicit part
    ! ghost data
    allocate( vbr1(4,ay,az),vbr2(4,ay,az) )
    allocate( vbs1(4,ay,az),vbs2(4,ay,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(ax-3:ax,:,:)
      vbs1 = v(1:4,:,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(ax-3:ax,:,:)
      vbr2 = v(1:4,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-3:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:4,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,az
    do j=1,ay
      dv(1,j,k) = sum(op%ar(1:4,1)*vbr1(1:4,j,k))+sum(op%ar(5:9,1)*v(1:5,j,k))
      dv(2,j,k) = sum(op%ar(1:3,2)*vbr1(2:4,j,k))+sum(op%ar(4:9,2)*v(1:6,j,k))
      dv(3,j,k) = sum(op%ar(1:2,3)*vbr1(3:4,j,k))+sum(op%ar(3:9,3)*v(1:7,j,k))
      dv(4,j,k) = sum(op%ar(1:1,4)*vbr1(4:4,j,k))+sum(op%ar(2:9,4)*v(1:8,j,k))
    end do
    end do
    do k=1,az
    do j=1,ay
      do i=5,ax-4
        dv(i,j,k) = sum(op%ar(1:9,i)*v(i-4:i+4,j,k))
      end do
    end do
    end do
    do k=1,az
    do j=1,ay
      dv(ax-3,j,k) = sum(op%ar(1:8,ax-3)*v(ax-7:ax,j,k))+sum(op%ar(9:9,ax-3)*vbr2(1:1,j,k))
      dv(ax-2,j,k) = sum(op%ar(1:7,ax-2)*v(ax-6:ax,j,k))+sum(op%ar(8:9,ax-2)*vbr2(1:2,j,k))
      dv(ax-1,j,k) = sum(op%ar(1:6,ax-1)*v(ax-5:ax,j,k))+sum(op%ar(7:9,ax-1)*vbr2(1:3,j,k))
      dv(ax  ,j,k) = sum(op%ar(1:5,ax  )*v(ax-4:ax,j,k))+sum(op%ar(6:9,ax  )*vbr2(1:4,j,k))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(ax,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ay,az), dvo(4,ay,az,0:np-1) )
      forall(i=1:2,j=1:ay,k=1:az) 
        dvop(i,j,k) = dv(i,j,k)
        dvop(i+2,j,k) = dv(ax+i-2,j,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ay, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ay, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:2)*dvo(3:4,j,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,3:4)*dvo(1:2,j,k,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ay, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1x_r4

  function eval_compact_op1y_r4(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=4, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r4), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in y'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
    ! explicit part
    ! ghost data
    allocate( vbr1(ax,4,az),vbr2(ax,4,az) )
    allocate( vbs1(ax,4,az),vbs2(ax,4,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,ay-3:ay,:)
      vbs1 = v(:,1:4,:)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,ay-3:ay,:)
      vbr2 = v(:,1:4,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-3:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:4,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do k=1,az
    do i=1,ax
      dv(i,1,k) = sum(op%ar(1:4,1)*vbr1(i,1:4,k))+sum(op%ar(5:9,1)*v(i,1:5,k))
      dv(i,2,k) = sum(op%ar(1:3,2)*vbr1(i,2:4,k))+sum(op%ar(4:9,2)*v(i,1:6,k))
      dv(i,3,k) = sum(op%ar(1:2,3)*vbr1(i,3:4,k))+sum(op%ar(3:9,3)*v(i,1:7,k))
      dv(i,4,k) = sum(op%ar(1:1,4)*vbr1(i,4:4,k))+sum(op%ar(2:9,4)*v(i,1:8,k))
    end do
    end do
    do k=1,az
    do j=5,ay-4
    do i=1,ax
        dv(i,j,k) = sum(op%ar(1:9,j)*v(i,j-4:j+4,k))
    end do
    end do
    end do
    do k=1,az
    do i=1,ax
      dv(i,ay-3,k) = sum(op%ar(1:8,ay-3)*v(i,ay-7:ay,k))+sum(op%ar(9:9,ay-3)*vbr2(i,1:1,k))
      dv(i,ay-2,k) = sum(op%ar(1:7,ay-2)*v(i,ay-6:ay,k))+sum(op%ar(8:9,ay-2)*vbr2(i,1:2,k))
      dv(i,ay-1,k) = sum(op%ar(1:6,ay-1)*v(i,ay-5:ay,k))+sum(op%ar(7:9,ay-1)*vbr2(i,1:3,k))
      dv(i,ay  ,k) = sum(op%ar(1:5,ay  )*v(i,ay-4:ay,k))+sum(op%ar(6:9,ay  )*vbr2(i,1:4,k))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,ay,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ax,az), dvo(4,ax,az,0:np-1) )
      forall(i=1:ax,j=1:2,k=1:az) 
        dvop(j,i,k) = dv(i,j,k)
        dvop(j+2,i,k) = dv(i,ay+j-2,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:2)*dvo(3:4,i,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,3:4)*dvo(1:2,i,k,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1y_r4

  function eval_compact_op1z_r4(op,v,vb1,vb2,dv1,dv2) result(dv) ! nor=4, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r4), intent(in) :: op
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(size(v,1),size(v,2),size(v,3)) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    integer :: nb,nsr,i,j,k
    integer :: ax,ay,az,nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v
!      print *,'null op in z'
      RETURN
    END IF
    ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
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
    np =  op%np
    ! explicit part
    ! ghost data
    allocate( vbr1(ax,ay,4),vbr2(ax,ay,4) )
    allocate( vbs1(ax,ay,4),vbs2(ax,ay,4) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,:,az-3:az)
      vbs1 = v(:,:,1:4)
      nsr = size(vbs1)
      CALL startCOMM()
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
      CALL endCOMM()
    else if( op%periodic ) then
      vbr1 = v(:,:,az-3:az)
      vbr2 = v(:,:,1:4)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-3:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:4)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    do j=1,ay
    do i=1,ax
      dv(i,j,1) = sum(op%ar(1:4,1)*vbr1(i,j,1:4))+sum(op%ar(5:9,1)*v(i,j,1:5))
      dv(i,j,2) = sum(op%ar(1:3,2)*vbr1(i,j,2:4))+sum(op%ar(4:9,2)*v(i,j,1:6))
      dv(i,j,3) = sum(op%ar(1:2,3)*vbr1(i,j,3:4))+sum(op%ar(3:9,3)*v(i,j,1:7))
      dv(i,j,4) = sum(op%ar(1:1,4)*vbr1(i,j,4:4))+sum(op%ar(2:9,4)*v(i,j,1:8))
    end do
    end do
    do k=5,az-4
    do j=1,ay
    do i=1,ax
        dv(i,j,k) = sum(op%ar(1:9,k)*v(i,j,k-4:k+4))
    end do
    end do
    end do
    do j=1,ay
    do i=1,ax
      dv(i,j,az-3) = sum(op%ar(1:8,az-3)*v(i,j,az-7:az))+sum(op%ar(9:9,az-3)*vbr2(i,j,1:1))
      dv(i,j,az-2) = sum(op%ar(1:7,az-2)*v(i,j,az-6:az))+sum(op%ar(8:9,az-2)*vbr2(i,j,1:2))
      dv(i,j,az-1) = sum(op%ar(1:6,az-1)*v(i,j,az-5:az))+sum(op%ar(7:9,az-1)*vbr2(i,j,1:3))
      dv(i,j,az  ) = sum(op%ar(1:5,az  )*v(i,j,az-4:az))+sum(op%ar(6:9,az  )*vbr2(i,j,1:4))
    end do
    end do
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,az)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) return
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) return
    ! parallel solver
      allocate( dvop(4,ax,ay), dvo(4,ax,ay,0:np-1) )
      forall(i=1:ax,j=1:ay,k=1:2) 
        dvop(k,i,j) = dv(i,j,k)
        dvop(k+2,i,j) = dv(i,j,az+k-2)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
         call startCOMM();call startCustom(1)
         call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
         call endCOMM();call endCustom(1)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, ay )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:2)*dvo(3:4,i,j,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,3:4)*dvo(1:2,i,j,op%hi))
     case( 2 ) ! mpi_gather/scatter
        call startCOMM();call startCustom(1)
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        call endCOMM();call endCustom(1)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, ay )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
      end select
      deallocate( dvop, dvo )
  end function eval_compact_op1z_r4

  subroutine setup_compact_ops(cops,patch,com)
    implicit none
    class(compact_type), intent(out) :: cops
    type(patch_type), intent(in) :: patch
    type(comm_type), intent(in) :: com
!    type(line_type) :: xline,yline,zline
    type(mesh1_type) :: xmsh,ymsh,zmsh
    type(comm1_type) :: xcom,ycom,zcom
    type(compact_weight), pointer :: weight
    logical(c_bool) :: null_op,null_opx,null_opy,null_opz
    logical :: spew = .FALSE.
    integer :: n,spec !,directcom
    ! get 1D from 3D types
    ! in terms of Miranda block and comm types
!     xmsh = mesh1_type(patch%coordsys,patch%ax,patch%nx,patch%dx,patch%x1,patch%xn,patch%bx1,patch%bxn)
     xmsh = mesh1_type(patch%ax,patch%nx,patch%dx,patch%x1,patch%xn,patch%bx1,patch%bxn)
     xcom = comm1_type(patch%periodicx,com%xcom,com%xcom_np,com%xcom_id, &
       com%xcom_lo,com%xcom_hi,com%xrange)
!     xline = line_type(patch%coordsys,patch%ax,patch%nx,patch%px,patch%dx, &
!       patch%x1,patch%xn,patch%bx1,patch%bxn,patch%isymx,patch%periodicx, &
!       com%xcom,com%xcom_np,com%xcom_id,com%xcom_lo,com%xcom_hi,com%xrange )
!     ymsh = mesh1_type(patch%coordsys,patch%ay,patch%ny,patch%dy,patch%y1,patch%yn,patch%by1,patch%byn)
     ymsh = mesh1_type(patch%ay,patch%ny,patch%dy,patch%y1,patch%yn,patch%by1,patch%byn)
     ycom = comm1_type(patch%periodicy,com%ycom,com%ycom_np,com%ycom_id, &
       com%ycom_lo,com%ycom_hi,com%yrange)
!     yline = line_type(patch%coordsys,patch%ay,patch%ny,patch%py,patch%dy, &
!       patch%y1,patch%yn,patch%by1,patch%byn,patch%isymy,patch%periodicy, &
!       com%ycom,com%ycom_np,com%ycom_id,com%ycom_lo,com%ycom_hi,com%yrange )
!     zmsh = mesh1_type(patch%coordsys,patch%az,patch%nz,patch%dz,patch%z1,patch%zn,patch%bz1,patch%bzn)
     zmsh = mesh1_type(patch%az,patch%nz,patch%dz,patch%z1,patch%zn,patch%bz1,patch%bzn)
     zcom = comm1_type(patch%periodicz,com%zcom,com%zcom_np,com%zcom_id, &
       com%zcom_lo,com%zcom_hi,com%zrange)
!     zline = line_type(patch%coordsys,patch%az,patch%nz,patch%pz,patch%dz, &
!       patch%z1,patch%zn,patch%bz1,patch%bzn,patch%isymz,patch%periodicz, &
!       com%zcom,com%zcom_np,com%zcom_id,com%zcom_lo,com%zcom_hi,com%zrange )
    ! unknown/periodic (default)
    cops%mbc = 0 ; cops%nop = 1
    ! symmetric-antisymmetric pairs
    if( xmsh%bc1=='SYMM' ) then ; cops%mbc(1,1,1:2) = [1,-1] ; cops%nop(1) = 2 ; endif
    if( xmsh%bcn=='SYMM' ) then ; cops%mbc(2,1,1:2) = [1,-1] ; cops%nop(1) = 2 ; endif
    if( ymsh%bc1=='SYMM' ) then ; cops%mbc(1,2,1:2) = [1,-1] ; cops%nop(2) = 2 ; endif
    if( ymsh%bcn=='SYMM' ) then ; cops%mbc(2,2,1:2) = [1,-1] ; cops%nop(2) = 2 ; endif
    if( zmsh%bc1=='SYMM' ) then ; cops%mbc(1,3,1:2) = [1,-1] ; cops%nop(3) = 2 ; endif
    if( zmsh%bcn=='SYMM' ) then ; cops%mbc(2,3,1:2) = [1,-1] ; cops%nop(3) = 2 ; endif
    ! extended rhs
    if( xmsh%bc1=='XTND' ) cops%mbc(1,1,1:2) = 2
    if( xmsh%bcn=='XTND' ) cops%mbc(2,1,1:2) = 2
    if( ymsh%bc1=='XTND' ) cops%mbc(1,2,1:2) = 2
    if( ymsh%bcn=='XTND' ) cops%mbc(2,2,1:2) = 2
    if( zmsh%bc1=='XTND' ) cops%mbc(1,3,1:2) = 2
    if( zmsh%bcn=='XTND' ) cops%mbc(2,3,1:2) = 2
    ! grid spacing
    cops%dx = xmsh%d ; cops%dy = ymsh%d ; cops%dz = zmsh%d
!    ! operation control (internal)
!    cops%control%null_opx = ( patch%nx < 4 .or. cops%control%null_opx )
!    cops%control%null_opy = ( patch%ny < 4 .or. cops%control%null_opy )
!    cops%control%null_opz = ( patch%nz < 4 .or. cops%control%null_opz )
    ! operation control (from input)
    cops%control%null_opx = ( patch%nx < 4 .or. zerodx )
    cops%control%null_opy = ( patch%ny < 4 .or. zerody )
    cops%control%null_opz = ( patch%nz < 4 .or. zerodz )
    cops%control%directcom = directcom
    ! local copies
!    directcom = cops%control%directcom
    null_opx = cops%control%null_opx ; null_opy = cops%control%null_opy ; null_opz = cops%control%null_opz
    ! make sure weights are set up
!    spew = (com%patcom_id == master) ! .false.
    if( .not. stencils_set ) call setup_stencils(spew)
    ! d1
    spec = cops%control%d1spec
    if( spec < 1 .or. spec > nderiv1 ) then  ! should be a fatal error
      null_op = .true.
      print *,'## invalid d1spec ##'
    else
      null_op = .false.
      weight => compact_weight_d1(spec)
    endif
    do n=1,cops%nop(1)
      cops%d1x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d1x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d1x(n)%directcom = directcom
      if( spew .and. .not. cops%d1x(n)%null_op ) print *,n,'setting up compact d1x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d1y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d1y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d1y(n)%directcom = directcom
      if( spew .and. .not. cops%d1y(n)%null_op ) print *,n,'setting up compact d1y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d1z(n)%null_op = null_op
      if( null_op ) cycle
     call cops%d1z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d1z(n)%directcom = directcom
      if( spew .and. .not. cops%d1z(n)%null_op ) print *,n,'setting up compact d1z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d2
    spec = cops%control%d2spec
    if( spec < 1 .or. spec > nderiv2 ) then
      null_op = .true.
      print *,'## invalid d2spec ##'
    else
      null_op = .false.
      weight => compact_weight_d2(spec)
    endif
    do n=1,cops%nop(1)
      cops%d2x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d2x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d2x(n)%directcom = directcom
      if( spew .and. .not. cops%d2x(n)%null_op ) print *,n,'setting up compact d2x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d2y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d2y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d2y(n)%directcom = directcom
      if( spew .and. .not. cops%d2y(n)%null_op ) print *,n,'setting up compact d2y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d2z(n)%null_op = null_op
      if( null_op ) cycle
     call cops%d2z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d2z(n)%directcom = directcom
      if( spew .and. .not. cops%d2z(n)%null_op ) print *,n,'setting up compact d2z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d4
    spec = cops%control%d4spec
    if( spec < 1 .or. spec > nderiv4 ) then
      null_op = .true.
      print *,'## invalid d4spec ##'
    else
      null_op = .false.
      weight => compact_weight_d4(spec)
    endif
    do n=1,cops%nop(1)
      cops%d4x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d4x(n)%directcom = directcom
      if( spew .and. .not. cops%d4x(n)%null_op ) print *,n,'setting up compact d4x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d4y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d4y(n)%directcom = directcom
      if( spew .and. .not. cops%d4y(n)%null_op ) print *,n,'setting up compact d4y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d4z(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d4z(n)%directcom = directcom
      if( spew .and. .not. cops%d4z(n)%null_op ) print *,n,'setting up compact d4z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d8
    spec = cops%control%d8spec
    if( spec < 1 .or. spec > nderiv8 ) then
      null_op = .true.
      print *,'## invalid d8spec ##'
    else
      null_op = .false.
      weight => compact_weight_d8(spec)
    endif
    do n=1,cops%nop(1)
      cops%d8x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d8x(n)%directcom = directcom
      if( spew .and. .not. cops%d8x(n)%null_op ) print *,n,'setting up compact d8x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d8y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d8y(n)%directcom = directcom
      if( spew .and. .not. cops%d8y(n)%null_op ) print *,n,'setting up compact d8y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d8z(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d8z(n)%directcom = directcom
      if( spew .and. .not. cops%d8z(n)%null_op ) print *,n,'setting up compact d8z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! gf
    spec = cops%control%gfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid gfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%gfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%gfx(n)%directcom = directcom
      if( spew .and. .not. cops%gfx(n)%null_op ) print *,n,'setting up compact gfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%gfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%gfy(n)%directcom = directcom
      if( spew .and. .not. cops%gfy(n)%null_op ) print *,n,'setting up compact gfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%gfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%gfz(n)%directcom = directcom
      if( spew .and. .not. cops%gfz(n)%null_op ) print *,n,'setting up compact gfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! sf
    spec = cops%control%sfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid sfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%sfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%sfx(n)%directcom = directcom
      if( spew .and. .not. cops%sfx(n)%null_op ) print *,n,'setting up compact sfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%sfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%sfy(n)%directcom = directcom
      if( spew .and. .not. cops%sfy(n)%null_op ) print *,n,'setting up compact sfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%sfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%sfz(n)%directcom = directcom
      if( spew .and. .not. cops%sfz(n)%null_op ) print *,n,'setting up compact sfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! tf
    spec = cops%control%tfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid tfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%tfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%tfx(n)%directcom = directcom
      if( spew .and. .not. cops%tfx(n)%null_op ) print *,n,'setting up compact tfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%tfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%tfy(n)%directcom = directcom
      if( spew .and. .not. cops%tfy(n)%null_op ) print *,n,'setting up compact tfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%tfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%tfz(n)%directcom = directcom
      if( spew .and. .not. cops%tfz(n)%null_op ) print *,n,'setting up compact tfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! imfc
    spec = cops%control%imfcspec
    if( spec < 1 .or. spec > nimfc ) then
      null_op = .true.
      print *,'## invalid imfcspec ##'
    else
      null_op = .false.
      weight => compact_weight_imfc(spec)
    endif
    do n=1,cops%nop(1)
      cops%imfcx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imfcx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%imfcx(n)%directcom = directcom
      if( spew .and. .not. cops%imfcx(n)%null_op ) print *,n,'setting up compact imfcx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%imfcy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imfcy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%imfcy(n)%directcom = directcom
      if( spew .and. .not. cops%imfcy(n)%null_op ) print *,n,'setting up compact imfcy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%imfcz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imfcz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%imfcz(n)%directcom = directcom
      if( spew .and. .not. cops%imfcz(n)%null_op ) print *,n,'setting up compact imfcz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! imcf
    spec = cops%control%imcfspec
    if( spec < 1 .or. spec > nimcf ) then
      null_op = .true.
      print *,'## invalid imcfspec ##'
    else
      null_op = .false.
      weight => compact_weight_imcf(spec)
    endif
    do n=1,cops%nop(1)
      cops%imcfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imcfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%imcfx(n)%directcom = directcom
      if( spew .and. .not. cops%imcfx(n)%null_op ) print *,n,'setting up compact imcfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%imcfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imcfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%imcfy(n)%directcom = directcom
      if( spew .and. .not. cops%imcfy(n)%null_op ) print *,n,'setting up compact imcfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%imcfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%imcfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%imcfz(n)%directcom = directcom
      if( spew .and. .not. cops%imcfz(n)%null_op ) print *,n,'setting up compact imcfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! isr  ! right shift
    spec = cops%control%isrspec
    if( spec < 1 .or. spec > nis ) then
      null_op = .true.
      print *,'## invalid isrspec ##'
    else
       null_op = .false.
       weight => compact_weight_ish(spec)
    endif
    do n=1,cops%nop(1)
      cops%isrx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%isrx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%isrx(n)%directcom = directcom
      if( spew .and. .not. cops%isrx(n)%null_op ) print *,n,'setting up compact isrx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%isry(n)%null_op = null_op
      if( null_op ) cycle
      call cops%isry(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%isry(n)%directcom = directcom
      if( spew .and. .not. cops%isry(n)%null_op ) print *,n,'setting up compact isry ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%isrz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%isrz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%isrz(n)%directcom = directcom
      if( spew .and. .not. cops%isrz(n)%null_op ) print *,n,'setting up compact isrz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! isl  ! left shift
    spec = cops%control%islspec
    if( spec < 1 .or. spec > nis ) then
      null_op = .true.
      print *,'## invalid islspec ##'
    else
       null_op = .false.
       weight => compact_weight_ish(spec)
    endif
    do n=1,cops%nop(1)
      cops%islx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%islx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%islx(n)%directcom = directcom
      if( spew .and. .not. cops%islx(n)%null_op ) print *,n,'setting up compact islx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%isly(n)%null_op = null_op
      if( null_op ) cycle
      call cops%isly(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%isly(n)%directcom = directcom
      if( spew .and. .not. cops%isly(n)%null_op ) print *,n,'setting up compact isly ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%islz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%islz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%islz(n)%directcom = directcom
      if( spew .and. .not. cops%islz(n)%null_op ) print *,n,'setting up compact islz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! amrcf
    spec = cops%control%amrcfspec
    if( spec < 1 .or. spec > nfamr ) then  ! should be a fatal error
      null_op = .true.
      print *,'## invalid amrcf ##'
    else
      null_op = .false.
      weight => compact_weight_famr(spec)
    endif
    do n=1,cops%nop(1)
      cops%amrcfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%amrcfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%amrcfx(n)%directcom = directcom
      if( spew .and. .not. cops%amrcfx(n)%null_op ) print *,n,'setting up compact amrcfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%amrcfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%amrcfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%amrcfy(n)%directcom = directcom
      if( spew .and. .not. cops%amrcfy(n)%null_op ) print *,n,'setting up compact amrcfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%amrcfz(n)%null_op = null_op
      if( null_op ) cycle
     call cops%amrcfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%amrcfz(n)%directcom = directcom
      if( spew .and. .not. cops%amrcfz(n)%null_op ) print *,n,'setting up compact amrcfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! amrfc
    spec = cops%control%amrfcspec
    if( spec < 1 .or. spec > nfamr ) then  ! should be a fatal error
      null_op = .true.
      print *,'## invalid amrfc ##'
    else
      null_op = .false.
      weight => compact_weight_famr(spec)
    endif
    do n=1,cops%nop(1)
      cops%amrfcx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%amrfcx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%amrfcx(n)%directcom = directcom
      if( spew .and. .not. cops%amrfcx(n)%null_op ) print *,n,'setting up compact amrfcx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%amrfcy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%amrfcy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%amrfcy(n)%directcom = directcom
      if( spew .and. .not. cops%amrfcy(n)%null_op ) print *,n,'setting up compact amrfcy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%amrfcz(n)%null_op = null_op
      if( null_op ) cycle
     call cops%amrfcz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%amrfcz(n)%directcom = directcom
      if( spew .and. .not. cops%amrfcz(n)%null_op ) print *,n,'setting up compact amrfcz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
  end subroutine setup_compact_ops

  subroutine remove_compact_ops(compact)
    class(compact_type), intent(out) :: compact
    continue
  end subroutine remove_compact_ops

end module LES_compact

