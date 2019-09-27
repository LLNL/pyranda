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
module LES_stencils
  USE iso_c_binding
  IMPLICIT NONE
  
  logical :: stencils_set=.false., verbose=.false.
  
  type compact_descriptor
    CHARACTER(KIND=c_char,LEN=8) :: name='void',complicity='illicit'
    integer(c_int) :: inorder=0,bcorder=0
  contains
    procedure :: info => write_compact_description
  end type compact_descriptor

  TYPE compact_weight
    INTEGER(c_int) :: nol,nor,nir,ncl,ncr,nci ! stencil
    INTEGER(c_int) :: nst  ! staggering, -1 = face to center, +1 = center to face
    INTEGER(c_int) :: nbc1,nbc2  ! telescoped boundary points
    INTEGER(c_int) :: null_option ! what to do if null_op: 0=zero, 1=copy, 2: = ??
    LOGICAL(c_bool) :: implicit_op  ! invert lhs; if nol == 0, implicit_op = .false.
    REAL(KIND=c_double) :: sh  ! shift from input grid
    type(compact_descriptor) :: description
    REAL(KIND=c_double), DIMENSION(:), ALLOCATABLE :: ali  ! lhs interior/periodic weights
    REAL(KIND=c_double), DIMENSION(:), ALLOCATABLE :: ari  ! rhs interior/periodic weights
    REAL(KIND=c_double), DIMENSION(:,:,:), ALLOCATABLE :: alb1,alb2  ! lhs telescoped boundary weights
    REAL(KIND=c_double), DIMENSION(:,:,:), ALLOCATABLE :: arb1,arb2  ! rhs telescoped boundary weights
  contains
    procedure :: info => write_weight_info
  END TYPE compact_weight
  
  INTEGER(c_int), PARAMETER :: nderiv1=2
  TYPE(compact_weight), target :: compact_weight_d1(nderiv1)
  INTEGER(c_int), PARAMETER :: nderiv2=1
  TYPE(compact_weight), target :: compact_weight_d2(nderiv2)
  INTEGER(c_int), PARAMETER :: nderiv4=1
  TYPE(compact_weight), target :: compact_weight_d4(nderiv4)
  INTEGER(c_int), PARAMETER :: nderiv8=1
  TYPE(compact_weight), target :: compact_weight_d8(nderiv8)
  INTEGER(c_int), PARAMETER :: nfilter=8
  TYPE(compact_weight), target :: compact_weight_ff(nfilter)
  INTEGER(c_int), PARAMETER :: nimcf=1
  TYPE(compact_weight), target :: compact_weight_imcf(nimcf)
  INTEGER(c_int), PARAMETER :: nimfc=1
  TYPE(compact_weight), target :: compact_weight_imfc(nimfc)
  INTEGER(c_int), PARAMETER :: nis=2 ! usually left and right shifts
  TYPE(compact_weight), target :: compact_weight_ish(nis)
  REAL(KIND=c_double), DIMENSION(nis) :: shift=[ -1.0_c_double/3.0_c_double, 1.0_c_double/3.0_c_double ]
  INTEGER(c_int), PARAMETER :: nfamr=2
  TYPE(compact_weight), target :: compact_weight_famr(nfamr)
  
  interface lower_symm_weights
    module procedure lower_symm_weights_int, lower_symm_weights_gen, lower_symm_weights_both
  end interface

  interface upper_symm_weights
    module procedure upper_symm_weights_int, upper_symm_weights_gen, upper_symm_weights_both
  end interface

contains

  subroutine write_compact_description(descriptor)
    class(compact_descriptor), intent(in) :: descriptor
    IF (.FALSE.) &
    print '(a,a9,a,a9,a,i3,a,i3)', &
             ' name =',descriptor%name, &
             ' , type =',descriptor%complicity, &
             ' , interior order =',descriptor%inorder, &
             ' , boundary order =',descriptor%bcorder
  end subroutine write_compact_description
  
  subroutine write_weight_info(weight)
    class(compact_weight), intent(in) :: weight
    character(len=10) :: lformat,rformat
    if( .not. stencils_set ) then
      print *, '## compact weights are not set ##'
      return
    endif
    call weight%description%info()
    write(lformat,'(a,i2,a)') '(',weight%ncl,'es14.6)'
    write(rformat,'(a,i2,a)') '(',weight%ncr,'es14.6)'
    print '(3(a,i3),a,f8.4)', &
             ' left stencil =',weight%ncl, &
             ' , right stencil =',weight%ncr, &
             ' , staggering =',weight%nst, &
             ' , shift =',weight%sh
    print *,'ali'
    print lformat,weight%ali
    print *,'alb1'
    print lformat,weight%alb1(:,:,:)
    print *,'alb2'
    print lformat,weight%alb2(:,:,:)
    print *,'ari'
    print rformat,weight%ari
    print *,'arb1'
    print rformat,weight%arb1(:,:,:)
    print *,'arb2'
    print rformat,weight%arb2(:,:,:)
  end subroutine write_weight_info
  
  subroutine setup_stencils(verb) ! (sh)
    implicit none
!    REAL(KIND=c_double), dimension(:), intent(in) :: sh
    logical, intent(in) :: verb
    integer :: s
!    verbose = verb
    if( verbose) print *,'available compact weights:'
    call c10d1(compact_weight_d1(1))  ! 1st derivative
    call c10d2(compact_weight_d2(1))  ! 2nd derivative
    call e4d4( compact_weight_d4(1))  ! 4th derivative (explicit)
    call c10d8(compact_weight_d8(1))  ! 8th derivative
    call c8ff9(compact_weight_ff(1))  ! 9/10 filter
    call c8ff8(compact_weight_ff(2))  ! 8/10 filter
    call c8ff7(compact_weight_ff(3))  ! 7/10 filter
    call c8ff6(compact_weight_ff(4))  ! 2/3 filter
    call cgft4(compact_weight_ff(5))  ! gaussian filter
    call cgfs4(compact_weight_ff(6))  ! gaussian filter w/ symmetry bcs
    call cgfc4(compact_weight_ff(7))  ! gaussian filter w/ constant bcs
    call ctfs4(compact_weight_ff(8))  ! tophat filter
    call c10imcf(compact_weight_imcf(1))  ! center-to-face midpoint interpolation
    call c10imfc(compact_weight_imfc(1))  ! face-to-center midpoint interpolation
!    do s = 1,size(sh)
!      if( s > nis ) exit
    do s = 1,nis
      call c10ish(compact_weight_ish(s),shift(s))  ! interpolation w/ arbitrary shift
    end do
    call cfamrcf(compact_weight_famr(1))  ! coarse-to-fine AMR filter
    call cfamrfc(compact_weight_famr(2))  ! fine-to-coarse AMR filter
!    print *,' '
    stencils_set = .true.
  end subroutine setup_stencils

  subroutine c10d1(d1)  ! setup compact 10th order 1st derivative weights
    implicit none
    TYPE(compact_weight), intent(out) :: d1
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 3, nir = 3
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    allocate( d1%ali(2*nol+1) )
    allocate( d1%ari(2*nor+1) )
    allocate( d1%alb1(2*nol+1,nbc1,-1:2), d1%alb2(2*nol+1,nbc2,-1:2) )
    allocate( d1%arb1(2*nor+1,nbc1,-1:2), d1%arb2(2*nor+1,nbc2,-1:2) )
    d1%description = compact_descriptor('deriv1st','explicit',10,3)
    if( implicit_op ) d1%description%complicity = 'implicit'
    if( verbose) call d1%description%info()
    d1%nol = nol  ! half-width of lhs stencil
    d1%nor = nor  ! half-width of outer rhs stencil
    d1%nir = nir  ! half-width of inner rhs stencil
    d1%ncl = ncl  ! total width of lhs stencil
    d1%ncr = ncr  ! total width of rhs stencil
    d1%nci = nci  ! parallel overlap of lhs stencil 
    d1%nst = nst  ! shift for staggered solution
    d1%nbc1 = nbc1  ! number of custom lower boundary points
    d1%nbc2 = nbc2  ! number of custom upper boundary points
    d1%null_option = 0 ! zero if null_op
    d1%implicit_op = implicit_op  ! lhs ignored if false
    if( nol == 0 ) d1%implicit_op = .false.
    d1%sh = 0.0d0  ! shift
  ! interior weights
    d1%ali = [ 0.45D0, 4.5D0, 9.0D0, 4.5D0, 0.45D0 ]
    d1%ari = [ -0.015D0, -1.515D0, -6.375D0, 0.0D0, 6.375D0, 1.515D0, 0.015D0 ]
  ! one-sided weights, band-diagonal style
    d1%alb1(:,1,0) = [ 0.0D0         , 0.0D0         , 4.725D0       , 9.45D0        , 0.0D0 ]
    d1%alb1(:,2,0) = [ 0.0D0         , 1.94578125D0  , 7.783125D0    , 1.94578125D0  , 0.0D0 ]
    d1%alb1(:,3,0) = [ 0.2964375D0   , 4.743D0       , 10.67175D0    , 4.743D0       , 0.2964375D0 ]
    d1%alb1(:,4,0) = [ 0.451390625D0 , 4.63271875D0  , 9.38146875D0  , 4.63271875D0  , 0.451390625D0 ]
    d1%alb2(:,1,0) = d1%alb1(ncl:1:-1,4,0)
    d1%alb2(:,2,0) = d1%alb1(ncl:1:-1,3,0)
    d1%alb2(:,3,0) = d1%alb1(ncl:1:-1,2,0)
    d1%alb2(:,4,0) = d1%alb1(ncl:1:-1,1,0)
    d1%arb1(:,1,0) = [ 0.0D0     , 0.0D0          , 0.0D0          , -11.8125D0  ,  9.45D0        , 2.3625D0      , 0.0D0   ]
    d1%arb1(:,2,0) = [ 0.0D0     , 0.0D0          , -5.83734375D0  , 0.0D0       ,  5.83734375D0  , 0.0D0         , 0.0D0   ]
    d1%arb1(:,3,0) = [ 0.0D0     , -1.23515625D0  , -7.905D0       , 0.0D0       ,  7.905D0       , 1.23515625D0  , 0.0D0   ]
    d1%arb1(:,4,0) = [ -0.015D0  , -1.53D0        , -6.66984375D0  , 0.0D0       ,  6.66984375D0  , 1.53D0        , 0.015D0 ]
    d1%arb2(:,1,0) = -d1%arb1(ncr:1:-1,4,0)
    d1%arb2(:,2,0) = -d1%arb1(ncr:1:-1,3,0)
    d1%arb2(:,3,0) = -d1%arb1(ncr:1:-1,2,0)
    d1%arb2(:,4,0) = -d1%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(d1%alb1(:,:,+1),d1%arb1(:,:,+1),d1%ali,d1%ari,nol,nor,0,0,-1,+1)
    call lower_symm_weights(d1%alb1(:,:,-1),d1%arb1(:,:,-1),d1%ali,d1%ari,nol,nor,0,0,+1,-1)
    call upper_symm_weights(d1%alb2(:,:,+1),d1%arb2(:,:,+1),d1%ali,d1%ari,nol,nor,0,0,-1,+1)
    call upper_symm_weights(d1%alb2(:,:,-1),d1%arb2(:,:,-1),d1%ali,d1%ari,nol,nor,0,0,+1,-1)
  ! with extended boundary data
    d1%alb1(:,1,2) = [ 0.0D0 , 0.0D0  , 6.0D0 , 0.0D0  , 0.0D0 ]
    d1%alb1(:,2,2) = [ 0.0D0 , 2.25D0 , 6.0D0 , 2.25D0 , 0.0D0 ]
    d1%alb1(:,3,2) = [ 0.3D0 , 3.0D0  , 6.0D0 , 3.0D0  , 0.3D0 ]
    d1%alb1(:,4,2) = [ 0.3D0 , 3.0D0  , 6.0D0 , 3.0D0  , 0.3D0 ]
    d1%alb2(:,1,2) = d1%alb1(ncl:1:-1,4,2)
    d1%alb2(:,2,2) = d1%alb1(ncl:1:-1,3,2)
    d1%alb2(:,3,2) = d1%alb1(ncl:1:-1,2,2)
    d1%alb2(:,4,2) = d1%alb1(ncl:1:-1,1,2)
    d1%arb1(:,1,2) = [ -0.1D0   , 0.9D0   , -4.5D0    , 0.0D0 , 4.5D0    , -0.9D0 , 0.1D0     ]
    d1%arb1(:,2,2) = [ 0.0125D0 , -0.3D0  , -4.6875D0 , 0.0D0 , 4.6875D0 , 0.3D0  , -0.0125D0 ]
    d1%arb1(:,3,2) = [ -0.01D0  , -1.01D0 , -4.25D0   , 0.0D0 , 4.25D0   , 1.01D0 , 0.01D0    ]
    d1%arb1(:,4,2) = [ -0.01D0  , -1.01D0 , -4.25D0   , 0.0D0 , 4.25D0   , 1.01D0 , 0.01D0    ]
    d1%arb2(:,1,2) = -d1%arb1(ncr:1:-1,4,2)
    d1%arb2(:,2,2) = -d1%arb1(ncr:1:-1,3,2)
    d1%arb2(:,3,2) = -d1%arb1(ncr:1:-1,2,2)
    d1%arb2(:,4,2) = -d1%arb1(ncr:1:-1,1,2)
  end subroutine c10d1
  
  subroutine c10d2(d2)  ! setup compact 10th order 2nd derivative weights
    implicit none
    TYPE(compact_weight), intent(out) :: d2
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 3, nir = 3
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    d2%description = compact_descriptor('deriv2nd','explicit',10,3)
    if( implicit_op ) d2%description%complicity = 'implicit'
    if( verbose) call d2%description%info()
    allocate( d2%ali(2*nol+1) )
    allocate( d2%ari(2*nor+1) )
    allocate( d2%alb1(2*nol+1,nbc1,-1:2), d2%alb2(2*nol+1,nbc2,-1:2) )
    allocate( d2%arb1(2*nor+1,nbc1,-1:2), d2%arb2(2*nor+1,nbc2,-1:2) )
    d2%nol = nol
    d2%nor = nor
    d2%nir = nir
    d2%ncl = ncl
    d2%ncr = ncr
    d2%nci = nci
    d2%nst = nst
    d2%nbc1 = nbc1
    d2%nbc2 = nbc2
    d2%null_option = 0 ! zero if null_op
    d2%implicit_op = implicit_op
    if( nol == 0 ) d2%implicit_op = .false.
    d2%sh = 0.0d0  ! shift
  ! interior weights
    d2%ali = [ 387.0D0 , 6012.0D0 , 16182.0D0 , 6012.0D0 , 387.0D0 ]
    d2%ari = [ 79.0D0  , 4671.0D0 , 9585.0D0  , -28670.0D0 , 9585.0D0 , 4671.0D0 , 79.0D0 ]
  ! one-sided weights, band-diagonal style
    d2%alb1(:,1,0) = [ 0.0D0   , 0.0D0    , 1.0D0     , 11.0D0   , 0.0D0   ]
    d2%alb1(:,2,0) = [ 0.0D0   , 1.0D0    , 10.0D0    , 1.0D0    , 0.0D0   ]
    d2%alb1(:,3,0) = [ 23.0D0  , 688.0D0  , 2358.0D0  , 688.0D0  , 23.0D0  ]
    d2%alb1(:,4,0) = [ 387.0D0 , 6012.0D0 , 16182.0D0 , 6012.0D0 , 387.0D0 ]
    d2%alb2(:,1,0) = d2%alb1(ncl:1:-1,4,0)
    d2%alb2(:,2,0) = d2%alb1(ncl:1:-1,3,0)
    d2%alb2(:,3,0) = d2%alb1(ncl:1:-1,2,0)
    d2%alb2(:,4,0) = d2%alb1(ncl:1:-1,1,0)
    d2%arb1(:,1,0) = [ 0.0D0   , 0.0D0    , 0.0D0     , 13.0D0     , -27.0D0  , 15.0D0   , -1.0D0  ]
    d2%arb1(:,2,0) = [ 0.0D0   , 0.0D0    , 12.0D0    , -24.0D0    , 12.0D0   , 0.0D0    , 0.0D0   ]
    d2%arb1(:,3,0) = [ 0.0D0   , 465.0D0  , 1920.0D0  , -4770.0D0  , 1920.0D0 , 465.0D0  , 0.0D0   ]
    d2%arb1(:,4,0) = [ 79.0D0  , 4671.0D0 , 9585.0D0  , -28670.0D0 , 9585.0D0 , 4671.0D0 , 79.0D0  ]
    d2%arb2(:,1,0) = d2%arb1(ncr:1:-1,4,0)
    d2%arb2(:,2,0) = d2%arb1(ncr:1:-1,3,0)
    d2%arb2(:,3,0) = d2%arb1(ncr:1:-1,2,0)
    d2%arb2(:,4,0) = d2%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(d2%alb1(:,:,+1),d2%arb1(:,:,+1),d2%ali,d2%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(d2%alb1(:,:,-1),d2%arb1(:,:,-1),d2%ali,d2%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(d2%alb2(:,:,+1),d2%arb2(:,:,+1),d2%ali,d2%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(d2%alb2(:,:,-1),d2%arb2(:,:,-1),d2%ali,d2%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    d2%alb1(:,1,2) = [ 0.0D0   , 0.0D0   , 18.0D0   , 0.0D0   , 0.0D0   ]
    d2%alb1(:,2,2) = [ 0.0D0   , 4.05D0  , 17.1D0   , 4.05D0  , 0.0D0   ]
    d2%alb1(:,3,2) = [ 0.387D0 , 6.012D0 , 16.182D0 , 6.012D0 , 0.387D0 ]
    d2%alb1(:,4,2) = [ 0.387D0 , 6.012D0 , 16.182D0 , 6.012D0 , 0.387D0 ]
    d2%alb2(:,1,2) = d2%alb1(ncl:1:-1,4,2)
    d2%alb2(:,2,2) = d2%alb1(ncl:1:-1,3,2)
    d2%alb2(:,3,2) = d2%alb1(ncl:1:-1,2,2)
    d2%alb2(:,4,2) = d2%alb1(ncl:1:-1,1,2)
    d2%arb1(:,1,2) = [ 0.2D0     , -2.7D0  , 27.0D0    , -49.0D0   , 27.0D0    , -2.7D0  , 0.2D0     ]
    d2%arb1(:,2,2) = [ -0.0575D0 , 2.295D0 , 16.5375D0 , -37.55D0  , 16.5375D0 , 2.295D0 , -0.0575D0 ]
    d2%arb1(:,3,2) = [ 0.079D0   , 4.671D0 , 9.585D0   , -28.670D0 , 9.585D0   , 4.671D0 , 0.079D0   ]
    d2%arb1(:,4,2) = [ 0.079D0   , 4.671D0 , 9.585D0   , -28.670D0 , 9.585D0   , 4.671D0 , 0.079D0   ]
    d2%arb2(:,1,2) = d2%arb1(ncr:1:-1,4,2)
    d2%arb2(:,2,2) = d2%arb1(ncr:1:-1,3,2)
    d2%arb2(:,3,2) = d2%arb1(ncr:1:-1,2,2)
    d2%arb2(:,4,2) = d2%arb1(ncr:1:-1,1,2)
  end subroutine c10d2

  subroutine c10d8(d8)  ! setup compact 10th (?) order 8th derivative weights
    implicit none
    TYPE(compact_weight), intent(out) :: d8
! Interior ---------------------------------------------------------------------
    REAL(c_double), PARAMETER :: zeta  =    29.0D0
    REAL(c_double), PARAMETER :: alpha =    14.0D0
    REAL(c_double), PARAMETER :: beta  =     1.5D0
    REAL(c_double), PARAMETER :: aa    =  4200.0D0
    REAL(c_double), PARAMETER :: bb    = -3360.0D0
    REAL(c_double), PARAMETER :: cc    =  1680.0D0
    REAL(c_double), PARAMETER :: dd    =  -480.0D0
    REAL(c_double), PARAMETER :: ee    =    60.0D0
! Symmetry boundaries ----------------------------------------------------------
    REAL(c_double), PARAMETER :: alpha2 =  beta + alpha
    
    REAL(c_double), PARAMETER :: zeta1  =  alpha + zeta
    REAL(c_double), PARAMETER :: alpha1 =  beta + alpha
           
    REAL(c_double), PARAMETER :: dd4   =  ee + dd
    
    REAL(c_double), PARAMETER :: cc3   =  dd + cc
    REAL(c_double), PARAMETER :: bb3   =  ee + bb
    
    REAL(c_double), PARAMETER :: bb2   =  cc + bb 
    REAL(c_double), PARAMETER :: aa2   =  dd + aa
    REAL(c_double), PARAMETER :: b22   =  ee + bb
    
    REAL(c_double), PARAMETER :: aa1   =  bb + aa
    REAL(c_double), PARAMETER :: bb1   =  cc + bb 
    REAL(c_double), PARAMETER :: cc1   =  dd + cc 
    REAL(c_double), PARAMETER :: dd1   =  ee + dd         
!------------------------------------------------------------------------------- 
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    d8%description = compact_descriptor('deriv8th','explicit',10,0)
    if( implicit_op ) d8%description%complicity = 'implicit'
    if( verbose) call d8%description%info()
    allocate( d8%ali(2*nol+1) )
    allocate( d8%ari(2*nor+1) )
    allocate( d8%alb1(2*nol+1,nbc1,-1:2), d8%alb2(2*nol+1,nbc2,-1:2) )
    allocate( d8%arb1(2*nor+1,nbc1,-1:2), d8%arb2(2*nor+1,nbc2,-1:2) )
    d8%nol = nol
    d8%nor = nor
    d8%nir = nir
    d8%ncl = ncl
    d8%ncr = ncr
    d8%nci = nci
    d8%nst = nst
    d8%nbc1 = nbc1
    d8%nbc2 = nbc2
    d8%null_option = 0 ! zero if null_op
    d8%implicit_op = implicit_op
    if( nol == 0 ) d8%implicit_op = .false.
    d8%sh = 0.0d0  ! shift
  ! interior weights
    d8%ali = [         beta, alpha, zeta, alpha, beta         ]
    d8%ari = [ ee, dd,   cc,    bb,   aa,    bb,   cc, dd, ee ]
  ! one-sided weights, band-diagonal style
    d8%alb1(:,1,0) = [ 0.0D0,  0.0D0,  zeta1, alpha1, beta ]
    d8%alb1(:,2,0) = [ 0.0D0, alpha2,   zeta,  alpha, beta ]
    d8%alb1(:,3,0) = [  beta,  alpha,   zeta,  alpha, beta ]
    d8%alb1(:,4,0) = [  beta,  alpha,   zeta,  alpha, beta ]
    d8%alb2(:,1,0) = d8%alb1(ncl:1:-1,4,0)
    d8%alb2(:,2,0) = d8%alb1(ncl:1:-1,3,0)
    d8%alb2(:,3,0) = d8%alb1(ncl:1:-1,2,0)
    d8%alb2(:,4,0) = d8%alb1(ncl:1:-1,1,0)
    d8%arb1(:,1,0) = [ 0.0D0, 0.0D0, 0.0D0,  0.0D0,  aa1,  bb1, cc1, dd1, ee ]
    d8%arb1(:,2,0) = [ 0.0D0, 0.0D0, 0.0D0,    bb2,  aa2,  b22,  cc,  dd, ee ]
    d8%arb1(:,3,0) = [ 0.0D0, 0.0D0,   cc3,    bb3,   aa,   bb,  cc,  dd, ee ]
    d8%arb1(:,4,0) = [ 0.0D0,   dd4,    cc,     bb,   aa,   bb,  cc,  dd, ee ]
    d8%arb2(:,1,0) = d8%arb1(ncr:1:-1,4,0)
    d8%arb2(:,2,0) = d8%arb1(ncr:1:-1,3,0)
    d8%arb2(:,3,0) = d8%arb1(ncr:1:-1,2,0)
    d8%arb2(:,4,0) = d8%arb1(ncr:1:-1,1,0)
! DIAGNOSTICS--------
!    WRITE(6,*) "Interior:", SUM(d8%ali        ),' - ',SUM(d8%ari        ),' = ',SUM(d8%ali)         - SUM(d8%ari)
!    WRITE(6,*) "      j4:", SUM(d8%alb1(:,4,0)),' - ',SUM(d8%arb1(:,4,0)),' = ',SUM(d8%alb1(:,4,0)) - SUM(d8%arb1(:,4,0))
!    WRITE(6,*) "      j3:", SUM(d8%alb1(:,3,0)),' - ',SUM(d8%arb1(:,3,0)),' = ',SUM(d8%alb1(:,3,0)) - SUM(d8%arb1(:,3,0))
!    WRITE(6,*) "      j2:", SUM(d8%alb1(:,2,0)),' - ',SUM(d8%arb1(:,2,0)),' = ',SUM(d8%alb1(:,2,0)) - SUM(d8%arb1(:,2,0))
!    WRITE(6,*) "      j1:", SUM(d8%alb1(:,1,0)),' - ',SUM(d8%arb1(:,1,0)),' = ',SUM(d8%alb1(:,1,0)) - SUM(d8%arb1(:,1,0))
!--------------------
  ! symmetric weights
    call lower_symm_weights(d8%alb1(:,:,+1),d8%arb1(:,:,+1),d8%ali,d8%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(d8%alb1(:,:,-1),d8%arb1(:,:,-1),d8%ali,d8%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(d8%alb2(:,:,+1),d8%arb2(:,:,+1),d8%ali,d8%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(d8%alb2(:,:,-1),d8%arb2(:,:,-1),d8%ali,d8%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    d8%alb1(:,1,2) = [ 0.0D0 ,  0.0D0  , 60.0D0 ,  0.0D0 , 0.0D0 ]
    d8%alb1(:,2,2) = [ 0.0D0 , 20.0D0  , 20.0D0 , 20.0D0 , 0.0D0 ]
    d8%alb1(:,3,2) = d8%ali
    d8%alb1(:,4,2) = d8%ali
    d8%alb2(:,1,2) = d8%alb1(ncl:1:-1,4,2)
    d8%alb2(:,2,2) = d8%alb1(ncl:1:-1,3,2)
    d8%alb2(:,3,2) = d8%alb1(ncl:1:-1,2,2)
    d8%alb2(:,4,2) = d8%alb1(ncl:1:-1,1,2)
    d8%arb1(:,1,2) = d8%ari
    d8%arb1(:,2,2) = d8%ari
    d8%arb1(:,3,2) = d8%ari
    d8%arb1(:,4,2) = d8%ari
    d8%arb2(:,1,2) = d8%arb1(ncr:1:-1,4,2)
    d8%arb2(:,2,2) = d8%arb1(ncr:1:-1,3,2)
    d8%arb2(:,3,2) = d8%arb1(ncr:1:-1,2,2)
    d8%arb2(:,4,2) = d8%arb1(ncr:1:-1,1,2)
  end subroutine c10d8
  
  subroutine c8ff9(ff)  ! setup compact 8th order 9/10 filter weights
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterc9','explicit',8,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = [ 1.6688D-1 , 6.6624D-1 , 1.0D0     , 6.6624D-1 , 1.6688D-1 ]
    ff%ari = [ -5.0D-6   , 4.0D-5    , 1.6674D-1 , 6.6652D-1 , 9.9965D-1 , 6.6652D-1 , 1.6674D-1 , 4.0D-5 , -5.0D-6 ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,2,0) = [ 0.0D0     , 4.997D-1  , 1.0D0 , 4.997D-1  , 0.0D0 ]
    ff%alb1(:,3,0) = [ 1.6688D-1 , 6.6624D-1 , 1.0D0 , 6.6624D-1 , 1.6688D-1 ]
    ff%alb1(:,4,0) = [ 1.6688D-1 , 6.6624D-1 , 1.0D0 , 6.6624D-1 , 1.6688D-1 ]
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0  , 0.0D0     , 0.0D0     , 1.0D0     , 0.0D0     , 0.0D0     , 0.0D0  , 0.0D0  ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0  , 0.0D0     , 4.9985D-1 , 9.997D-1  , 4.9985D-1 , 0.0D0     , 0.0D0  , 0.0D0  ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0  , 1.668D-1  , 6.6656D-1 , 9.9952D-1 , 6.6656D-1 , 1.668D-1  , 0.0D0  , 0.0D0  ]
    ff%arb1(:,4,0) = [ 0.0D0 , 4.0D-5 , 1.6672D-1 , 6.6652D-1 , 9.9968D-1 , 6.6652D-1 , 1.6672D-1 , 4.0D-5 , 0.0D0  ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    ff%alb1(:,1,2) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,2,2) = [ 0.0D0 , 0.1950165953348364D0 , 0.6099668093303272D0 , 0.1950165953348364D0 , 0.0D0 ]
!    ff%alb1(:,2,2) = [ 0.0D0     , 4.997D-1  , 1.0D0 , 4.997D-1  , 0.0D0 ]
    ff%alb1(:,3,2) = ff%ali
    ff%alb1(:,4,2) = ff%ali
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
!    ff%arb1(:,1,2) = [  0.01171875D0 , -0.03125D0 , -0.046875D0 , 0.28125D0 , 0.5703125D0 , &  ! 0.5 + T"=0
!                        0.28125D0 , -0.046875D0 , -0.03125D0 ,  0.01171875D0]
    ff%arb1(:,1,2) = [ -0.00390625D0 ,  0.03125D0 , -0.109375D0 , 0.21875D0 , 0.7265625D0 , &  ! 0.7
                        0.21875D0 , -0.109375D0 ,  0.03125D0 , -0.00390625D0]
    ff%arb1(:,2,2) = [ -0.0008591156978932D0, &  ! 0.8
                        0.0068729255831458D0, &
                       -0.0240552395410097D0, &
                        0.2431270744168542D0, &
                        0.5498287104778058D0, &
                        0.2431270744168542D0, &
                       -0.0240552395410097D0, &
                        0.0068729255831458D0, &
                       -0.0008591156978932D0  ]
!    ff%arb1(:,1,2) = [ 0.0D0 , 0.0D0  , 0.0D0     , 0.0D0     , 1.0D0     , 0.0D0     , 0.0D0     , 0.0D0  , 0.0D0  ]
!    ff%arb1(:,2,2) = [ 0.0D0 , 0.0D0  , 0.0D0     , 4.9985D-1 , 9.997D-1  , 4.9985D-1 , 0.0D0     , 0.0D0  , 0.0D0  ]
    ff%arb1(:,3,2) = ff%ari
    ff%arb1(:,4,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine c8ff9

  subroutine c8ff8(ff)  ! setup compact 8th order 8/10 filter weights
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterc8','explicit',8,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = [ 1.7136D-1 , 6.5728D-1 , 1.0D0     , 6.5728D-1 , 1.7136D-1 ]
    ff%ari = [ -1.1D-4   , 8.8D-4    , 1.6828D-1 , 6.6344D-1 , 9.923D-1  , 6.6344D-1 , 1.6828D-1 , 8.8D-4 , -1.1D-4 ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0    ]
    ff%alb1(:,2,0) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0    ]
    ff%alb1(:,3,0) = [ 1.6666666666666667D-1, 6.5585333333333333D-1, 1.0D0 , 6.5585333333333333D-1, 1.6666666666666667D-1 ]
    ff%alb1(:,4,0) = [ 1.6872D-1 , 6.564D-1  , 1.0D0 , 6.564D-1 , 1.6872D-1 ]
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0      , 0.0D0      , 9.375D-1  , 2.5D-1     , -3.75D-1   , 2.5D-1  , -6.25D-2  ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0      , 6.25D-2    , 7.5D-1    , 3.75D-1    , -2.5D-1    , 6.25D-2 , 0.0D0     ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , 1.65315D-1 , 6.6126D-1  , 9.9189D-1 , 6.6126D-1  , 1.65315D-1 , 0.0D0   , 0.0D0     ]
    ff%arb1(:,4,0) = [ 0.0D0 , 3.85D-4 , 1.6641D-1  , 6.62175D-1 , 9.923D-1  , 6.62175D-1 , 1.6641D-1  , 3.85D-4 , 0.0D0     ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    ff%alb1(:,1,2) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0    ]
    ff%alb1(:,2,2) = [ 0.0D0     , 0.0D0     , 1.0D0 , 0.0D0     , 0.0D0    ]
    ff%alb1(:,3,2) = [ 1.6666666666666667D-1, 6.5585333333333333D-1, 1.0D0 , 6.5585333333333333D-1, 1.6666666666666667D-1 ]
    ff%alb1(:,4,2) = [ 1.6872D-1 , 6.564D-1  , 1.0D0 , 6.564D-1 , 1.6872D-1 ]
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = [ 0.0D0 , 0.0D0   , 0.0D0      , 0.0D0      , 9.375D-1  , 2.5D-1     , -3.75D-1   , 2.5D-1  , -6.25D-2  ]
    ff%arb1(:,2,2) = [ 0.0D0 , 0.0D0   , 0.0D0      , 6.25D-2    , 7.5D-1    , 3.75D-1    , -2.5D-1    , 6.25D-2 , 0.0D0     ]
    ff%arb1(:,3,2) = [ 0.0D0 , 0.0D0   , 1.65315D-1 , 6.6126D-1  , 9.9189D-1 , 6.6126D-1  , 1.65315D-1 , 0.0D0   , 0.0D0     ]
    ff%arb1(:,4,2) = [ 0.0D0 , 3.85D-4 , 1.6641D-1  , 6.62175D-1 , 9.923D-1  , 6.62175D-1 , 1.6641D-1  , 3.85D-4 , 0.0D0     ]
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine c8ff8

  subroutine c8ff7(ff)  ! setup compact 8th order 7/10 filter weights
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterc7','explicit',8,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = [ 1.95D-1  , 6.1D-1   , 1.0D0   , 6.1D-1   , 1.95D-1 ]
    ff%ari = [ -6.640625D-4 , 5.3125D-3 , 1.7640625D-1 , 6.471875D-1 , 9.53515625D-1 , &
                6.471875D-1 , 1.7640625D-1 , 5.3125D-3 , -6.640625D-4 ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,2,0) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,3,0) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,4,0) = [ 0.3D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.3D0 ]
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0     , 0.0D0    , 0.9375D0 ,  0.25D0 , -0.375D0  , 0.25D0   , -0.0625D0 ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0     , 0.0625D0 ,  0.75D0  , 0.375D0 , -0.25D0   , 0.0625D0 , 0.0D0     ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , -0.0625D0 , 0.25D0   , 0.625D0  , 0.25D0  , -0.0625D0 , 0.0D0    , 0.0D0     ]
    ff%arb1(:,4,0) = [ 0.0D0 , 0.025D0 , 0.15D0    , 0.375D0  , 0.5D0    , 0.375D0 , 0.15D0    , 0.025D0  , 0.0D0     ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    ff%alb1(:,1,2) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,2,2) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,3,2) = [ 0.0D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.0D0 ]
    ff%alb1(:,4,2) = [ 0.3D0     , 0.0D0    , 1.0D0 , 0.0D0     , 0.3D0 ]
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = [ 0.0D0 , 0.0D0   , 0.0D0     , 0.0D0    , 0.9375D0 ,  0.25D0 , -0.375D0  , 0.25D0   , -0.0625D0 ]
    ff%arb1(:,2,2) = [ 0.0D0 , 0.0D0   , 0.0D0     , 0.0625D0 ,  0.75D0  , 0.375D0 , -0.25D0   , 0.0625D0 , 0.0D0     ]
    ff%arb1(:,3,2) = [ 0.0D0 , 0.0D0   , -0.0625D0 , 0.25D0   , 0.625D0  , 0.25D0  , -0.0625D0 , 0.0D0    , 0.0D0     ]
    ff%arb1(:,4,2) = [ 0.0D0 , 0.025D0 , 0.15D0    , 0.375D0  , 0.5D0    , 0.375D0 , 0.15D0    , 0.025D0  , 0.0D0     ]
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine c8ff7

  subroutine c8ff6(ff)  ! setup compact 8th order 2/3 filter weights
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterc6','explicit',8,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = [ 2.2384D-1 , 5.5232D-1 , 1.0D0     , 5.5232D-1 , 2.2384D-1 ]
    ff%ari = [ -1.34D-3  , 1.072D-2  , 1.8632D-1 , 6.2736D-1 , 9.062D-1 , 6.2736D-1 , 1.8632D-1 , 1.072D-2 , -1.34D-3 ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,2,0) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,3,0) = [ 0.25D0   , 0.836D0 , 1.5D0 , 0.836D0 , 0.25D0   ]
    ff%alb1(:,4,0) = [ 1.912D-1 , 5.44D-1 , 1.0D0 , 5.44D-1 , 1.912D-1 ]
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0D0    , 0.9375D0 ,  0.25D0 , -0.375D0 , 0.25D0   , -0.0625D0 ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0625D0 ,  0.75D0  , 0.375D0 , -0.25D0  , 0.0625D0 , 0.0D0     ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , 0.2295D0 , 0.918D0  , 1.377D0  , 0.918D0 , 0.2295D0 , 0.0D0    , 0.0D0     ]
    ff%arb1(:,4,0) = [ 0.0D0 , 4.6D-3  , 1.636D-1 , 6.13D-1  , 9.08D-1  , 6.13D-1 , 1.636D-1 , 4.6D-3   , 0.0D0     ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    ff%alb1(:,1,2) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,2,2) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,3,2) = [ 0.25D0   , 0.836D0 , 1.5D0 , 0.836D0 , 0.25D0   ]
    ff%alb1(:,4,2) = [ 1.912D-1 , 5.44D-1 , 1.0D0 , 5.44D-1 , 1.912D-1 ]
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0D0    , 0.9375D0 ,  0.25D0 , -0.375D0 , 0.25D0   , -0.0625D0 ]
    ff%arb1(:,2,2) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0625D0 ,  0.75D0  , 0.375D0 , -0.25D0  , 0.0625D0 , 0.0D0     ]
    ff%arb1(:,3,2) = [ 0.0D0 , 0.0D0   , 0.2295D0 , 0.918D0  , 1.377D0  , 0.918D0 , 0.2295D0 , 0.0D0    , 0.0D0     ]
    ff%arb1(:,4,2) = [ 0.0D0 , 4.6D-3  , 1.636D-1 , 6.13D-1  , 9.08D-1  , 6.13D-1 , 1.636D-1 , 4.6D-3   , 0.0D0     ]
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine c8ff6

  subroutine cgft4(ff)  ! setup 4-point Gaussian with telescoped boundaries
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 0, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .false.
    DOUBLE PRECISION, PARAMETER :: agau = 3565.0D0/10368.0D0
    DOUBLE PRECISION, PARAMETER :: bgau = 3091.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: cgau = 1997.0D0/25920.0D0
    DOUBLE PRECISION, PARAMETER :: dgau = 149.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: egau = 107.0D0/103680.0D0
    DOUBLE PRECISION, PARAMETER :: a4gau = 17.0D0/48.0D0
    DOUBLE PRECISION, PARAMETER :: b4gau = 15.0D0/64.0D0
    DOUBLE PRECISION, PARAMETER :: c4gau = 7.0D0/96.0D0
    DOUBLE PRECISION, PARAMETER :: d4gau = 1.0D0/64.0D0
    DOUBLE PRECISION, PARAMETER :: a3gau = 31.0D0/64.0D0
    DOUBLE PRECISION, PARAMETER :: b3gau = 7.0D0/32.0D0
    DOUBLE PRECISION, PARAMETER :: c3gau = 5.0D0/128.0D0
    DOUBLE PRECISION, PARAMETER :: a2gau = 2.0D0/3.0D0
    DOUBLE PRECISION, PARAMETER :: b2gau = 1.0D0/6.0D0
    DOUBLE PRECISION, PARAMETER :: a1gau = 5.0D0/6.0D0
    DOUBLE PRECISION, PARAMETER :: b1gau = 1.0D0/6.0D0
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filtergt','explicit',0,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = 1.0D0
    ff%ari = [ egau ,  dgau , cgau , bgau , agau , bgau , cgau , dgau , egau ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = 1.0D0
    ff%alb1(:,2,0) = 1.0D0
    ff%alb1(:,3,0) = 1.0D0
    ff%alb1(:,4,0) = 1.0D0
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0D0    , a1gau    , b1gau  , 0.0D0     , 0.0D0   , 0.0D0 ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , b2gau    , a2gau    , b2gau  , 0.0D0     , 0.0D0   , 0.0D0 ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , c3gau    , b3gau    , a3gau    , b3gau  , c3gau     , 0.0D0   , 0.0D0 ]
    ff%arb1(:,4,0) = [ 0.0D0 , d4gau   , c4gau    , b4gau    , a4gau    , b4gau  , c4gau     , d4gau   , 0.0D0 ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    ff%alb1(:,1,2) = 1.0D0
    ff%alb1(:,2,2) = 1.0D0
    ff%alb1(:,3,2) = 1.0D0
    ff%alb1(:,4,2) = 1.0D0
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = ff%ari
    ff%arb1(:,2,2) = ff%ari
    ff%arb1(:,3,2) = ff%ari
    ff%arb1(:,4,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine cgft4

  subroutine c4ff3(ff)  ! setup compact 4th order 1/3 filter weights for AMR coarsening
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterc3','explicit',4,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = [ 0.22382150680239417D0, -0.55235698639521165D0, 1.0D0, -0.55235698639521165D0, 0.22382150680239417D0 ]
    ff%ari = [ 1.33956656568111D-3, 1.071653252544891D-2, 3.750786383907118D-2, 7.501572767814235D-2, 9.376965959767794D-2, &
      7.501572767814235D-2, 3.750786383907118D-2, 1.071653252544891D-2, 1.33956656568111D-3 ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,2,0) = [ 0.0D0    , 0.5D0   , 1.0D0 , 0.5D0   , 0.0D0    ]
    ff%alb1(:,3,0) = [ 0.1666666666666667D0, -0.5404638224773227D0, 1.0D0, -0.5404638224773227D0, 0.1666666666666667D0 ]
    ff%alb1(:,4,0) = [ 0.1925352813232545D0, -0.5373235933837273D0, 1.0D0, -0.5373235933837273D0, 0.1925352813232545D0 ]
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0D0    , 1.0D0 ,  0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.5D0    , 1.0D0  , 0.5D0 , 0.0D0 , 0.0D0 , 0.0D0 ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , 1.57753555236680D-2, 6.31014220946720D-2, 9.46521331420080D-2, &
      6.31014220946720D-2, 1.57753555236680D-2, 0.0D0, 0.0D0 ]
    ff%arb1(:,4,0) = [ 0.0D0, 4.8503652481102D-3, 2.91021914886613D-2, 7.27554787216534D-2, 9.70073049622046D-2, &
      7.27554787216534D-2, 2.91021914886613D-2, 4.8503652481102D-3, 0.0D0 ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    ff%alb1(:,1,2) = [ 0.0D0    , 0.0D0   , 1.0D0 , 0.0D0   , 0.0D0    ]
    ff%alb1(:,2,2) = [ 0.0D0    , 0.5D0   , 1.0D0 , 0.5D0   , 0.0D0    ]
    ff%alb1(:,3,2) = [ 0.1666666666666667D0, -0.5404638224773227D0, 1.0D0, -0.5404638224773227D0, 0.1666666666666667D0 ]
    ff%alb1(:,4,2) = [ 0.1925352813232545D0, -0.5373235933837273D0, 1.0D0, -0.5373235933837273D0, 0.1925352813232545D0 ]
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.0D0    , 1.0D0 ,  0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ]
    ff%arb1(:,2,2) = [ 0.0D0 , 0.0D0   , 0.0D0    , 0.5D0    , 1.0D0  , 0.5D0 , 0.0D0 , 0.0D0 , 0.0D0 ]
    ff%arb1(:,3,2) = [ 0.0D0 , 0.0D0   , 1.57753555236680D-2, 6.31014220946720D-2, 9.46521331420080D-2, &
      6.31014220946720D-2, 1.57753555236680D-2, 0.0D0, 0.0D0 ]
    ff%arb1(:,4,2) = [ 0.0D0, 4.8503652481102D-3, 2.91021914886613D-2, 7.27554787216534D-2, 9.70073049622046D-2, &
      7.27554787216534D-2, 2.91021914886613D-2, 4.8503652481102D-3, 0.0D0 ]
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine c4ff3

  subroutine cgfs4(ff)  ! setup 4-point Gaussian with symmetry boundaries
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 0, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .false.
    DOUBLE PRECISION, PARAMETER :: agau = 3565.0D0/10368.0D0
    DOUBLE PRECISION, PARAMETER :: bgau = 3091.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: cgau = 1997.0D0/25920.0D0
    DOUBLE PRECISION, PARAMETER :: dgau = 149.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: egau = 107.0D0/103680.0D0
    DOUBLE PRECISION, PARAMETER :: dg14 = dgau + egau
    DOUBLE PRECISION, PARAMETER :: cg13 = cgau + dgau
    DOUBLE PRECISION, PARAMETER :: bg23 = bgau + egau
    DOUBLE PRECISION, PARAMETER :: bg12 = bgau + cgau
    DOUBLE PRECISION, PARAMETER :: ag22 = agau + dgau
    DOUBLE PRECISION, PARAMETER :: bg32 = bgau + egau
    DOUBLE PRECISION, PARAMETER :: ag11 = agau + bgau
    DOUBLE PRECISION, PARAMETER :: bg21 = bgau + cgau
    DOUBLE PRECISION, PARAMETER :: cg31 = cgau + dgau
    DOUBLE PRECISION, PARAMETER :: dg41 = dgau + egau
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filtergs','explicit',0,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = 1.0D0
    ff%ari = [ egau ,  dgau , cgau , bgau , agau , bgau , cgau , dgau , egau ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = 1.0D0
    ff%alb1(:,2,0) = 1.0D0
    ff%alb1(:,3,0) = 1.0D0
    ff%alb1(:,4,0) = 1.0D0
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , 0.0D0   , ag11    , bg21     , cg31     , dg41     , egau ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , bg12    , ag22    , bg32     , cgau     , dgau     , egau ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , cg13    , bg23    , agau    , bgau     , cgau     , dgau     , egau ]
    ff%arb1(:,4,0) = [ 0.0D0 , dg14    , cgau    , bgau    , agau    , bgau     , cgau     , dgau     , egau ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    ff%alb1(:,1,2) = 1.0D0
    ff%alb1(:,2,2) = 1.0D0
    ff%alb1(:,3,2) = 1.0D0
    ff%alb1(:,4,2) = 1.0D0
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = ff%ari
    ff%arb1(:,2,2) = ff%ari
    ff%arb1(:,3,2) = ff%ari
    ff%arb1(:,4,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine cgfs4

  subroutine cgfc4(ff)  ! setup 4-point Gaussian with constant boundaries
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 0, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .false.
    DOUBLE PRECISION, PARAMETER :: agau = 3565.0D0/10368.0D0
    DOUBLE PRECISION, PARAMETER :: bgau = 3091.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: cgau = 1997.0D0/25920.0D0
    DOUBLE PRECISION, PARAMETER :: dgau = 149.0D0/12960.0D0
    DOUBLE PRECISION, PARAMETER :: egau = 107.0D0/103680.0D0
    DOUBLE PRECISION, PARAMETER :: dbdry = dgau + egau
    DOUBLE PRECISION, PARAMETER :: cbdry = cgau + dbdry
    DOUBLE PRECISION, PARAMETER :: bbdry = bgau + cbdry
    DOUBLE PRECISION, PARAMETER :: abdry = agau + bbdry
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filtergc','explicit',0,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = 1.0D0
    ff%ari = [ egau ,  dgau , cgau , bgau , agau , bgau , cgau , dgau , egau ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = 1.0D0
    ff%alb1(:,2,0) = 1.0D0
    ff%alb1(:,3,0) = 1.0D0
    ff%alb1(:,4,0) = 1.0D0
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , 0.0D0   , abdry    , bgau     , cgau     , dgau     , egau ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , bbdry   , agau     , bgau     , cgau     , dgau     , egau ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , cbdry   , bgau    , agau     , bgau     , cgau     , dgau     , egau ]
    ff%arb1(:,4,0) = [ 0.0D0 , dbdry   , cgau    , bgau     , agau    , bgau     , cgau     , dgau     , egau ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    ff%alb1(:,1,2) = 1.0D0
    ff%alb1(:,2,2) = 1.0D0
    ff%alb1(:,3,2) = 1.0D0
    ff%alb1(:,4,2) = 1.0D0
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = ff%ari
    ff%arb1(:,2,2) = ff%ari
    ff%arb1(:,3,2) = ff%ari
    ff%arb1(:,4,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine cgfc4

  subroutine ctfs4(ff)  ! setup 4-point tophat with symmetry boundaries
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 0, nor = 4, nir = 4
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .false.
    DOUBLE PRECISION, PARAMETER :: ahat = 1.0D0/9.0D0
    DOUBLE PRECISION, PARAMETER :: bhat = 1.0D0/9.0D0
    DOUBLE PRECISION, PARAMETER :: chat = 1.0D0/9.0D0
    DOUBLE PRECISION, PARAMETER :: dhat = 1.0D0/9.0D0
    DOUBLE PRECISION, PARAMETER :: ehat = 1.0D0/9.0D0
    DOUBLE PRECISION, PARAMETER :: dh14 = dhat + ehat
    DOUBLE PRECISION, PARAMETER :: ch13 = chat + dhat
    DOUBLE PRECISION, PARAMETER :: bh23 = bhat + ehat
    DOUBLE PRECISION, PARAMETER :: bh12 = bhat + chat
    DOUBLE PRECISION, PARAMETER :: ah22 = ahat + dhat
    DOUBLE PRECISION, PARAMETER :: bh32 = bhat + ehat
    DOUBLE PRECISION, PARAMETER :: ah11 = ahat + bhat
    DOUBLE PRECISION, PARAMETER :: bh21 = bhat + chat
    DOUBLE PRECISION, PARAMETER :: ch31 = chat + dhat
    DOUBLE PRECISION, PARAMETER :: dh41 = dhat + ehat
    INTEGER :: i,j,k
    ff%description = compact_descriptor('filterts','explicit',0,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    ff%ali = 1.0D0
    ff%ari = [ ehat , dhat , chat , bhat , ahat , bhat , chat , dhat , ehat ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = 1.0D0
    ff%alb1(:,2,0) = 1.0D0
    ff%alb1(:,3,0) = 1.0D0
    ff%alb1(:,4,0) = 1.0D0
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,4,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,4,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , 0.0D0   , ah11    , bh21     , ch31     , dh41     , ehat ]
    ff%arb1(:,2,0) = [ 0.0D0 , 0.0D0   , 0.0D0   , bh12    , ah22    , bh32     , chat     , dhat     , ehat ]
    ff%arb1(:,3,0) = [ 0.0D0 , 0.0D0   , ch13    , bh23    , ahat    , bhat     , chat     , dhat     , ehat ]
    ff%arb1(:,4,0) = [ 0.0D0 , dh14    , chat    , bhat    , ahat    , bhat     , chat     , dhat     , ehat ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,4,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,4,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    ff%alb1(:,1,2) = 1.0D0
    ff%alb1(:,2,2) = 1.0D0
    ff%alb1(:,3,2) = 1.0D0
    ff%alb1(:,4,2) = 1.0D0
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,4,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,4,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = ff%ari
    ff%arb1(:,2,2) = ff%ari
    ff%arb1(:,3,2) = ff%ari
    ff%arb1(:,4,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,4,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,4,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine ctfs4

  subroutine c10imcf(im)  ! setup compact 10th order center-to-face interpolation weights
    implicit none
    TYPE(compact_weight), intent(out) :: im
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 3, nir = 2
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 1, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    allocate( im%ali(ncl) )
    allocate( im%ari(ncr) )
    allocate( im%alb1(ncl,nbc1,-1:2), im%alb2(ncl,nbc2,-1:2) )
    allocate( im%arb1(ncr,nbc1,-1:2), im%arb2(ncr,nbc2,-1:2) )
    im%description = compact_descriptor('interpcf','explicit',10,2)
    if( implicit_op ) im%description%complicity = 'implicit'
    if( verbose) call im%description%info()
    im%nol = nol ; im%nor = nor ; im%nir = nir ; im%ncl = ncl ; im%ncr = ncr ; im%nci = nci
    im%nst = nst ; im%nbc1 = nbc1 ; im%nbc2 = nbc2
    im%null_option = 2 ! ?? if null_op
    im%implicit_op = implicit_op
    if( nol == 0 ) im%implicit_op = .false.
    im%sh = -0.5d0  ! shift
  ! interior weights
    im%ali = [ 0.01953125d0  , 0.234375d0    , 0.4921875d0   , 0.234375d0   , 0.01953125d0 ]
    im%ari = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0  , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
  ! one-sided weights, band-diagonal style
    im%alb1(:,1,0) = [ 0.0D0         , 0.0D0      , 1.0D0       , 0.0D0      , 0.0D0 ] ! explicit
    im%alb1(:,2,0) = [ 0.0D0         , 0.0d0      , 1.0d0       , 0.0d0      , 0.0D0 ] ! explicit
    im%alb1(:,3,0) = [ 0.0D0         , 0.1875d0   , 0.625d0     , 0.1875d0   , 0.0D0 ]
!    im%alb1(:,2,0) = [ 0.0D0         , 0.2d0      , 1.0d0       , 0.0d0      , 0.0D0 ]
!    im%alb1(:,3,0) = [ 0.01D0        , 0.27D0     , 0.63D0      , 0.21D0     , 0.0D0 ]
    im%alb1(:,4,0) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb2(:,1,0) = im%alb1(ncl:1:-1,4,0)
    im%alb2(:,2,0) = im%alb1(ncl:1:-1,3,0)
    im%alb2(:,3,0) = im%alb1(ncl:1:-1,2,0)
    im%alb2(:,4,0) = im%alb1(ncl:1:-1,1,0)
    im%arb1(:,1,0) = [ 0.0D0         , 0.0D0         , 0.0D0       , 1.5d0      , -0.5d0      , 0.0d0      ] ! extrapolate
    im%arb1(:,2,0) = [ 0.0D0         , 0.0D0         , 0.5D0       , 0.5D0      ,  0.0D0      , 0.0D0      ]
    im%arb1(:,3,0) = [ 0.0D0         , 0.03125d0     , 0.46875d0   , 0.46875d0  , 0.03125d0   , 0.0D0      ]
!    im%arb1(:,1,0) = [ 0.0D0         , 0.0D0         , 0.0D0        , 1.875d0      , -1.25d0       , 0.375d0       ] ! extrapolate
!    im%arb1(:,2,0) = [ 0.0D0         , 0.0D0         , 0.75D0       , 0.5D0        ,  -0.05D0      , 0.0D0         ]
!    im%arb1(:,3,0) = [ 0.0D0         , 0.07875d0     , 0.525d0      , 0.4725d0     ,  0.045d0      , -0.00125D0    ]
    im%arb1(:,4,0) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb2(:,1,0) = im%arb1(ncr:1:-1,4,0)
    im%arb2(:,2,0) = im%arb1(ncr:1:-1,3,0)
    im%arb2(:,3,0) = im%arb1(ncr:1:-1,2,0)
    im%arb2(:,4,0) = im%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(im%alb1(:,:,+1),im%arb1(:,:,+1),im%ali,im%ari,nol,nor,1,0,+1,+1)
    call lower_symm_weights(im%alb1(:,:,-1),im%arb1(:,:,-1),im%ali,im%ari,nol,nor,1,0,-1,-1)
    call upper_symm_weights(im%alb2(:,:,+1),im%arb2(:,:,+1),im%ali,im%ari,nol,nor,1,0,+1,+1)
    call upper_symm_weights(im%alb2(:,:,-1),im%arb2(:,:,-1),im%ali,im%ari,nol,nor,1,0,-1,-1)
  ! with extended boundary data -- needs work
    im%alb1(:,1,2) = [ 0.0D0         , 0.0D0      , 1.0D0       , 0.0D0    , 0.0D0 ]
    im%alb1(:,2,2) = [ 0.0D0         , 0.25d0     , 0.7d0       , 0.25d0   , 0.0D0 ]
    im%alb1(:,3,2) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb1(:,4,2) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb2(:,1,2) = im%alb1(ncl:1:-1,4,2)
    im%alb2(:,2,2) = im%alb1(ncl:1:-1,3,2)
    im%alb2(:,3,2) = im%alb1(ncl:1:-1,2,2)
    im%alb2(:,4,2) = im%alb1(ncl:1:-1,1,2)
    im%arb1(:,1,2) = [ 0.01171875D0  , -0.09765625D0 , 0.5859375D0  , 0.5859375d0  , -0.09765625d0 , 0.01171875d0  ]
    im%arb1(:,2,2) = [ -0.0015625D0  , 0.0546875d0   , 0.546875d0   , 0.546875d0   , 0.0546875d0   , -0.0015625D0  ]
    im%arb1(:,3,2) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb1(:,4,2) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb2(:,1,2) = im%arb1(ncr:1:-1,4,2)
    im%arb2(:,2,2) = im%arb1(ncr:1:-1,3,2)
    im%arb2(:,3,2) = im%arb1(ncr:1:-1,2,2)
    im%arb2(:,4,2) = im%arb1(ncr:1:-1,1,2)
  end subroutine c10imcf

  subroutine c10imfc(im)  ! setup compact 10th order center-to-face interpolation weights
    implicit none
    TYPE(compact_weight), intent(out) :: im
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 2, nir = 3
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = -1, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    INTEGER :: i,j,k
    im%description = compact_descriptor('interpfc','explicit',10,2)
    if( implicit_op ) im%description%complicity = 'implicit'
    if( verbose) call im%description%info()
    allocate( im%ali(ncl) )
    allocate( im%ari(ncr) )
    allocate( im%alb1(ncl,nbc1,-1:2), im%alb2(ncl,nbc2,-1:2) )
    allocate( im%arb1(ncr,nbc1,-1:2), im%arb2(ncr,nbc2,-1:2) )
    im%nol = nol ; im%nor = nor ; im%nir = nir ; im%ncl = ncl ; im%ncr = ncr ; im%nci = nci
    im%nst = nst ; im%nbc1 = nbc1 ; im%nbc2 = nbc2
    im%null_option = 3 ! ?? if null_op
    im%implicit_op = implicit_op
    if( nol == 0 ) im%implicit_op = .false.
    im%sh = 0.5d0  ! shift
  ! interior weights
    im%ali = [ 0.01953125d0  , 0.234375d0    , 0.4921875d0   , 0.234375d0   , 0.01953125d0 ]
    im%ari = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0  , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
  ! one-sided weights, band-diagonal style
    im%alb1(:,1,0) = [ 0.0D0         , 0.0D0      , 1.0D0       , 0.0D0      , 0.0D0 ]
    im%alb1(:,2,0) = [ 0.0D0         , 0.1875d0   , 0.625d0     , 0.1875d0   , 0.0D0 ]
    im%alb1(:,3,0) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb1(:,4,0) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb2(:,1,0) = im%alb1(ncl:1:-1,4,0)
    im%alb2(:,2,0) = im%alb1(ncl:1:-1,3,0)
    im%alb2(:,3,0) = im%alb1(ncl:1:-1,2,0)
    im%alb2(:,4,0) = im%alb1(ncl:1:-1,1,0)
    im%arb1(:,1,0) = [ 0.0D0         , 0.0D0         , 0.5D0        , 0.5d0      , 0.0d0       , 0.0d0       ]
    im%arb1(:,2,0) = [ 0.0D0         , 0.03125d0     , 0.46875d0    , 0.46875d0  , 0.03125d0   , 0.0D0       ]
    im%arb1(:,3,0) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb1(:,4,0) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb2(:,1,0) = im%arb1(ncr:1:-1,4,0)
    im%arb2(:,2,0) = im%arb1(ncr:1:-1,3,0)
    im%arb2(:,3,0) = im%arb1(ncr:1:-1,2,0)
    im%arb2(:,4,0) = im%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(im%alb1(:,:,+1),im%arb1(:,:,+1),im%ali,im%ari,nol,nor,0,1,+1,+1)
    call lower_symm_weights(im%alb1(:,:,-1),im%arb1(:,:,-1),im%ali,im%ari,nol,nor,0,1,-1,-1)
    call upper_symm_weights(im%alb2(:,:,+1),im%arb2(:,:,+1),im%ali,im%ari,nol,nor,0,1,+1,+1)
    call upper_symm_weights(im%alb2(:,:,-1),im%arb2(:,:,-1),im%ali,im%ari,nol,nor,0,1,-1,-1)
  ! with extended boundary data -- needs work
    im%alb1(:,1,2) = [ 0.0D0         , 0.0D0      , 1.0D0       , 0.0D0    , 0.0D0 ]
    im%alb1(:,2,2) = [ 0.0D0         , 0.25d0     , 0.7d0       , 0.25d0   , 0.0D0 ]
    im%alb1(:,3,2) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb1(:,4,2) = [ 0.01953125d0  , 0.234375d0 , 0.4921875d0 , 0.234375d0 , 0.01953125d0 ]
    im%alb2(:,1,2) = im%alb1(ncl:1:-1,4,2)
    im%alb2(:,2,2) = im%alb1(ncl:1:-1,3,2)
    im%alb2(:,3,2) = im%alb1(ncl:1:-1,2,2)
    im%alb2(:,4,2) = im%alb1(ncl:1:-1,1,2)
    im%arb1(:,1,2) = [ 0.01171875D0  , -0.09765625D0 , 0.5859375D0  , 0.5859375d0  , -0.09765625d0 , 0.01171875d0  ]
    im%arb1(:,2,2) = [ -0.0015625D0  , 0.0546875d0   , 0.546875d0   , 0.546875d0   , 0.0546875d0   , -0.0015625D0  ]
    im%arb1(:,3,2) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb1(:,4,2) = [ 0.001953125d0 , 0.087890625d0 , 0.41015625d0 , 0.41015625d0 , 0.087890625d0 , 0.001953125d0 ]
    im%arb2(:,1,2) = im%arb1(ncr:1:-1,4,2)
    im%arb2(:,2,2) = im%arb1(ncr:1:-1,3,2)
    im%arb2(:,3,2) = im%arb1(ncr:1:-1,2,2)
    im%arb2(:,4,2) = im%arb1(ncr:1:-1,1,2)
  end subroutine c10imfc

  subroutine c10ish(is,sin)  ! setup compact 10th order right-shift interpolation weights
    implicit none
    TYPE(compact_weight), intent(out) :: is
    real(kind=c_double), intent(in) :: sin  ! shift
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 3, nir = 3
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 4, nbc2 = 4
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    ! local working arrays
    real(kind=c_double), dimension(ncl) :: ali6,ali7,ali6b
    real(kind=c_double), dimension(ncr) :: ari6,ari7,ari6b
    real(kind=c_double), dimension(ncl,nbc1,-1:2) :: alb1
    real(kind=c_double), dimension(ncr,nbc1,-1:2) :: arb1
    real(kind=c_double), dimension(ncl,nbc2,-1:2) :: alb2
    real(kind=c_double), dimension(ncr,nbc2,-1:2) :: arb2
    real(kind=c_double) :: s
    INTEGER :: i,j,k
    is%description = compact_descriptor('interpsh','explicit',10,3)
    if( sin < 0.0 ) is%description%name = 'interpsl'
    if( sin > 0.0 ) is%description%name = 'interpsr'
    if( implicit_op ) is%description%complicity = 'implicit'
    if( verbose) call is%description%info()
    is%nol = nol ; is%nor = nor ; is%nir = nir ; is%ncl = ncl ; is%ncr = ncr ; is%nci = nci
    is%nst = nst ; is%nbc1 = nbc1 ; is%nbc2 = nbc2
    is%null_option = 4 ! ?? if null_op
    is%implicit_op = implicit_op
    if( nol == 0 ) is%implicit_op = .false.
    is%sh = sin  ! shift
    s = abs(sin) ! calculate for s > 0
  ! interior weights, 5 x 7 stencil
    ali7(1) = (1.0d0-s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)/1008.0d0
    ali7(2) = (3.0d0+s)*(4.0d0+s)*(5.0d0-3.0d0*s+s*s)/252.0d0
    ali7(3) = (16.0d0-s*s)*(5.0d0+s*s)/168.0d0
    ali7(4) = (3.0d0-s)*(4.0d0-s)*(5.0d0+3.0d0*s+s*s)/252.0d0
    ali7(5) = (1.0d0+s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)/1008.0d0
    ari7(1) = s*(1.0d0-s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)*(s-0.5d0)/15120.0d0
    ari7(2) = (1.0d0-s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)*(2.5d0-3.0d0*s-s*s)/2520.0d0
    ari7(3) = (4.0d0+s)*(1.0d0-s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)*(2.5d0+s)/1008.0d0
    ari7(4) = (9.0d0-s*s)*(16.0d0-s*s)*(2.5d0-s*s)/756.0d0
    ari7(5) = (4.0d0-s)*(1.0d0+s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)*(2.5d0-s)/1008.0d0
    ari7(6) = (1.0d0+s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)*(2.5d0+3.0d0*s-s*s)/2520.0d0
    ari7(7) = s*(1.0d0+s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)*(0.5d0+s)/15120.0d0
  ! interior weights, 5 x 6 stencil
    ali6(1) = (1.0d0+s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)/3024.0d0
    ali6(2) = (5.0d0-s)*(2.0d0+s)*(3.0d0+s)*(4.0d0+s)/756.0d0
    ali6(3) = (5.0d0-s)*(4.0d0-s)*(3.0d0+s)*(4.0d0+s)/504.0d0
    ali6(4) = (5.0d0-s)*(4.0d0-s)*(3.0d0-s)*(4.0d0+s)/756.0d0
    ali6(5) = (5.0d0-s)*(4.0d0-s)*(3.0d0-s)*(2.0d0-s)/3024.0d0
    ari6(1) = 0.0d0
    ari6(2) = (1.0d0-s)*ali6(5)/5.0d0
    ari6(3) = (4.0d0+s)*ali6(5)
    ari6(4) = (3.0d0+s)*ali6(4)/2.0d0
    ari6(5) = (2.0d0+s)*ali6(3)/3.0d0
    ari6(6) = (1.0d0+s)*ali6(2)/4.0d0
    ari6(7) = s*ali6(1)/5.0d0
  ! boundary weights, 5 x 6 stencil w/ extrapolation
    ali6b(5) = (1.0d0-s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)/3024.0d0
    ali6b(4) = (5.0d0+s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)/756.0d0
    ali6b(3) = (5.0d0+s)*(4.0d0+s)*(3.0d0-s)*(4.0d0-s)/504.0d0
    ali6b(2) = (5.0d0+s)*(4.0d0+s)*(3.0d0+s)*(4.0d0-s)/756.0d0
    ali6b(1) = (5.0d0+s)*(4.0d0+s)*(3.0d0+s)*(2.0d0+s)/3024.0d0
    ari6b(7) = 0.0d0
    ari6b(6) = (1.0d0+s)*ali6b(1)/5.0d0
    ari6b(5) = (4.0d0-s)*ali6b(1)
    ari6b(4) = (3.0d0-s)*ali6b(2)/2.0d0
    ari6b(3) = (2.0d0-s)*ali6b(3)/3.0d0
    ari6b(2) = (1.0d0-s)*ali6b(4)/4.0d0
    ari6b(1) = -s*ali6b(5)/5.0d0
  ! lower one-sided weights, band-diagonal style
    alb1 = 0.0d0 ; alb2 = 0.0d0 ; arb1 = 0.0d0 ; arb2 = 0.0d0
!    alb1(3,1,0) = 1.0d0 ! 1 x 2
!    alb1(3,1,0) = 1.0d0 ! 1 x 3
    alb1(3:4,1,0) = [ 1.0d0+s, 2.0d0-s ]/3.0d0 ! 2 x 3
    alb1(2,2,0) = 0.1d0+s*( 0.15d0+0.05d0*s)   ! 3 x 4
    alb1(3,2,0) = 0.6d0+s*( 0.10d0-0.10d0*s)
    alb1(4,2,0) = 0.3d0+s*(-0.25d0+0.05d0*s)
    alb1(:,3,0) = ali6 ! 5 x 6
    alb1(:,4,0) = ali7 ! 5 x 7
!    arb1(4:5,1,0) = [ 1.0d0-s , s ]
!    arb1(4:6,1,0) = [ (1.0d0-s)*(1.0d0-0.5d0*s), s*(2.0d0-s), -0.5d0*s*(1.0d0-s) ]
    arb1(4:6,1,0) = [ 2.0d0-s*(3.0d0-s), 4.0d0+2.0d0*s*(1.0d0-s), s*(1.0d0+s) ]/6.0d0
    arb1(3,2,0) = ( 6.0d0+s*(-11.0d0+s*( 6.0d0-s)))/60.0d0
    arb1(4,2,0) = (36.0d0+s*(-12.0d0+s*(-9.0d0+s*3.0d0)))/60.0d0
    arb1(5,2,0) = (18.0d0+s*( 21.0d0+s*( 0.0d0-s*3.0d0)))/60.0d0
    arb1(6,2,0) = (       s*(  2.0d0+s*( 3.0d0+s)))/60.0d0
    arb1(:,3,0) = ari6
    arb1(:,4,0) = ari7
  ! upper one-sided weights w/ extrapolation, band-diagonal style
!    alb2(3,4,0) = 1.0d0 ! 1 x 2
!    alb2(3,4,0) = 1.0d0 ! 1 x 3
    alb2(2:3,4,0) = [ 2.0d0+s, 1.0d0-s ]/3.0d0 ! 2 x 3
    alb2(4,3,0) = 0.1d0-s*( 0.15d0-0.05d0*s)   ! 3 x 4
    alb2(3,3,0) = 0.6d0-s*( 0.10d0+0.10d0*s)
    alb2(2,3,0) = 0.3d0-s*(-0.25d0-0.05d0*s)
    alb2(:,2,0) = ali6b ! 5 x 6
    alb2(:,1,0) = ali7  ! 5 x 7
!    arb2(3:4,4,0) = [ -s, 1.0d0+s ]  ! extrapolate on right boundary
!    arb2(2:4,4,0) = [ 0.5d0*s*(1.0d0+s), -s*(2.0d0+s), (1.0d0+s)*(1.0d0+0.5d0*s) ]
    arb2(2:4,4,0) = [ -s*(1.0d0-s), 4.0d0-2.0d0*s*(1.0d0+s), 2.0d0+s*(3.0d0+s) ]/6.0d0
    arb2(5,3,0) = ( 6.0d0-s*(-11.0d0-s*( 6.0d0+s)))/60.0d0
    arb2(4,3,0) = (36.0d0-s*(-12.0d0-s*(-9.0d0-s*3.0d0)))/60.0d0
    arb2(3,3,0) = (18.0d0-s*( 21.0d0-s*( 0.0d0+s*3.0d0)))/60.0d0
    arb2(2,3,0) = (      -s*(  2.0d0-s*( 3.0d0-s)))/60.0d0
    arb2(:,2,0) = ari6b
    arb2(:,1,0) = ari7
!   ! symmetric weights ! in general only rhs has discernible symmetry
!     call lower_symm_weights(alb1(:,:,+1),arb1(:,:,+1),ali7,ari7,nol,nor,0,0,+1,+1)
!     call lower_symm_weights(alb1(:,:,-1),arb1(:,:,-1),ali7,ari7,nol,nor,0,0,-1,-1)
!     call upper_symm_weights(alb2(:,:,+1),arb2(:,:,+1),ali7,ari7,nol,nor,0,0,+1,+1)
!     call upper_symm_weights(alb2(:,:,-1),arb2(:,:,-1),ali7,ari7,nol,nor,0,0,-1,-1)
  ! with extended boundary data -- needs work
    alb1(3,1,2) = 1.0D0
!     alb1(2,2,2) = (2.0d0+s)*(3.0d0+s)/42.0d0
!     alb1(3,2,2) = (4.0d0-s)*(3.0d0+s)/21.0d0
!     alb1(4,2,2) = (4.0d0-s)*(3.0d0-s)/42.0d0
    alb1(2,2,2) = ( 9.0d0-s*s)/42.0d0
    alb1(3,2,2) = (12.0d0+s*s)/21.0d0
    alb1(4,2,2) = ( 9.0d0-s*s)/42.0d0
    alb1(:,3,2) = ali7
    alb1(:,4,2) = ali7
    alb2(:,1,2) = alb1(:,4,2)
    alb2(:,2,2) = alb1(:,3,2)
    alb2(:,3,2) = alb1(:,2,2)
    alb2(:,4,2) = alb1(:,1,2)
!     arb1(2,1,2) = s*(3.0d0-s)*(2.0d0-s)*(1.0d0-s*s)/120.0d0
!     arb1(3,1,2) =-s*(3.0d0-s)*(1.0d0-s)*(4.0d0-s*s)/24.0d0
!     arb1(4,1,2) = (3.0d0-s)*(1.0d0-s*s)*(4.0d0-s*s)/12.0d0
!     arb1(5,1,2) = s*(3.0d0-s)*(1.0d0+s)*(4.0d0-s*s)/12.0d0
!     arb1(6,1,2) =-s*(3.0d0-s)*(2.0d0+s)*(1.0d0-s*s)/24.0d0
!     arb1(7,1,2) =         s*(4.0d0-s*s)*(1.0d0-s*s)/120.0d0
    arb1(1,1,2) =-s*(3.0d0-s)*(4.0d0-s*s)*(1.0d0-s*s)/720.0d0
    arb1(2,1,2) = s*(2.0d0-s)*(9.0d0-s*s)*(1.0d0-s*s)/120.0d0
    arb1(3,1,2) =-s*(1.0d0-s)*(9.0d0-s*s)*(4.0d0-s*s)/48.0d0
    arb1(4,1,2) = (1.0d0-s*s)*(9.0d0-s*s)*(4.0d0-s*s)/36.0d0
    arb1(5,1,2) = s*(1.0d0+s)*(9.0d0-s*s)*(4.0d0-s*s)/48.0d0
    arb1(6,1,2) =-s*(2.0d0+s)*(9.0d0-s*s)*(1.0d0-s*s)/120.0d0
    arb1(7,1,2) = s*(3.0d0+s)*(4.0d0-s*s)*(1.0d0-s*s)/720.0d0
!     arb1(2,2,2) = -s*(1.0d0-s)*(2.0d0-s)*(3.0d0-s)*(4.0d0-s)/2520.0d0
!     arb1(3,2,2) = (1.0d0-s)*(2.0d0-s)*(9.0d0-s*s)*(4.0d0-s)/504.0d0
!     arb1(4,2,2) = (4.0d0-s*s)*(9.0d0-s*s)*(4.0d0-s)/252.0d0
!     arb1(5,2,2) = (1.0d0+s)*(2.0d0+s)*(9.0d0-s*s)*(4.0d0-s)/252.0d0
!     arb1(6,2,2) = s*(1.0d0+s)*(2.0d0+s)*(3.0d0+s)*(4.0d0-s)/504.0d0
!     arb1(7,2,2) = -s*(1.0d0-s*s)*(2.0d0+s)*(3.0d0+s)/2520.0d0
    arb1(1,2,2) = s*(3.0d0-s)*(2.0d0-s)*(1.0d0-s*s)*(0.5d0-s)/2520.0d0
    arb1(2,2,2) =-s*(1.0d0-s)*(2.0d0-s)*(9.0d0-s*s)*(2.0d0-3.0d0*s)/1260.0d0
    arb1(3,2,2) =   (1.0d0-s)*(2.0d0-s)*(9.0d0-s*s)*(12.0d0-7.0d0*s-6.0d0*s*s)/1008.0d0
    arb1(4,2,2) =   (9.0d0-s*s)*(4.0d0-s*s)*(2.0d0-s*s)/126.0d0
    arb1(5,2,2) =   (1.0d0+s)*(2.0d0+s)*(9.0d0-s*s)*(12.0d0+7.0d0*s-6.0d0*s*s)/1008.0d0
    arb1(6,2,2) = s*(1.0d0+s)*(2.0d0+s)*(9.0d0-s*s)*(2.0d0+3.0d0*s)/1260.0d0
    arb1(7,2,2) =-s*(3.0d0+s)*(2.0d0+s)*(1.0d0-s*s)*(0.5d0+s)/2520.0d0
    arb1(:,3,2) = ari7
    arb1(:,4,2) = ari7
    arb2(:,1,2) = arb1(:,4,2)
    arb2(:,2,2) = arb1(:,3,2)
    arb2(:,3,2) = arb1(:,2,2)
    arb2(:,4,2) = arb1(:,1,2)
    ! symmetric weights ! impose only on rhs using extended weights
    arb1(:,:,1) = arb1(:,:,2)             ; arb1(:,:,-1) = arb1(:,:,2)
    arb2(:,:,1) = arb2(:,:,2)             ; arb2(:,:,-1) = arb2(:,:,2)
    call lower_symm_weights(alb1,arb1,nol,nor,0,0,+1,+1)
    call upper_symm_weights(alb2,arb2,nol,nor,0,0,+1,+1)
    alb1(:,:,1) = alb1(:,:,2)             ; alb1(:,:,-1) = alb1(:,:,2)
    alb2(:,:,1) = alb2(:,:,2)             ; alb2(:,:,-1) = alb2(:,:,2)
    allocate( is%ali(ncl) )
    allocate( is%ari(ncr) )
    allocate( is%alb1(ncl,nbc1,-1:2), is%alb2(ncl,nbc2,-1:2) )
    allocate( is%arb1(ncr,nbc1,-1:2), is%arb2(ncr,nbc2,-1:2) )
    if( sin >= 0.0d0 ) then
      is%ali = ali7 ; is%alb1 = alb1 ; is%alb2 = alb2
      is%ari = ari7 ; is%arb1 = arb1 ; is%arb2 = arb2
    else  ! reverse order
      is%ali = ali7(ncl:1:-1) ; is%ari = ari7(ncr:1:-1)
      is%alb1(:,1,:) = alb2(ncl:1:-1,4,:) ; is%alb2(:,1,:) = alb1(ncl:1:-1,4,:)
      is%alb1(:,2,:) = alb2(ncl:1:-1,3,:) ; is%alb2(:,2,:) = alb1(ncl:1:-1,3,:)
      is%alb1(:,3,:) = alb2(ncl:1:-1,2,:) ; is%alb2(:,3,:) = alb1(ncl:1:-1,2,:)
      is%alb1(:,4,:) = alb2(ncl:1:-1,1,:) ; is%alb2(:,4,:) = alb1(ncl:1:-1,1,:)
      is%arb1(:,1,:) = arb2(ncr:1:-1,4,:) ; is%arb2(:,1,:) = arb1(ncr:1:-1,4,:)
      is%arb1(:,2,:) = arb2(ncr:1:-1,3,:) ; is%arb2(:,2,:) = arb1(ncr:1:-1,3,:)
      is%arb1(:,3,:) = arb2(ncr:1:-1,2,:) ; is%arb2(:,3,:) = arb1(ncr:1:-1,2,:)
      is%arb1(:,4,:) = arb2(ncr:1:-1,1,:) ; is%arb2(:,4,:) = arb1(ncr:1:-1,1,:)
    endif
  end subroutine c10ish

!===================================================================================================
! INVERTIBLE AMR FILTER FOR COARSE-TO-FINE AND FINE-TO-COARSE OPERATIONS (fbar <---> f):
! alpha*fbar(i-1) + fbar(i) + alpha*fbar(i+1) = c*f(i-2) + b*f(i-1) + a*f(i) + b*f(i+1) + c*f(i+2)
! alpha = -0.0321826755129339
!     a =  0.4451523642186118
!     b =  0.2207614172195584
!     c =  0.0244797251582018
! The above coefficients match the transfer function of a gaussian filter of width 3*dx.
!===================================================================================================

  subroutine cfamrcf(cf)  ! setup compact coarse-to-fine for 3x refinement
    implicit none
    TYPE(compact_weight), intent(out) :: cf
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 2, nir = 2
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 3, nbc2 = 3  ! ?
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    REAL(c_double), PARAMETER :: two=2.0_c_double
    INTEGER :: i,j,k
    allocate( cf%ali(2*nol+1) )
    allocate( cf%ari(2*nor+1) )
    allocate( cf%alb1(2*nol+1,nbc1,-1:2), cf%alb2(2*nol+1,nbc2,-1:2) )
    allocate( cf%arb1(2*nor+1,nbc1,-1:2), cf%arb2(2*nor+1,nbc2,-1:2) )
    cf%description = compact_descriptor('famrcf13','explicit',0,0)
    if( implicit_op ) cf%description%complicity = 'implicit'
    call cf%description%info()
    cf%nol = nol  ! half-width of lhs stencil
    cf%nor = nor  ! half-width of outer rhs stencil
    cf%nir = nir  ! half-width of inner rhs stencil
    cf%ncl = ncl  ! total width of lhs stencil
    cf%ncr = ncr  ! total width of rhs stencil
    cf%nci = nci  ! parallel overlap of lhs stencil 
    cf%nst = nst  ! shift for staggered solution
    cf%nbc1 = nbc1  ! number of custom lower boundary points
    cf%nbc2 = nbc2  ! number of custom upper boundary points
    cf%null_option = 1 ! copy if null_op
    cf%implicit_op = implicit_op  ! lhs ignored if false
    if( nol == 0 ) cf%implicit_op = .false.
    cf%sh = 0.0d0  ! shift
  ! interior weights
    cf%ali = [ 0.0244797251582018D0, 0.2207614172195584D0, 0.4451523642186118D0, 0.2207614172195584D0, 0.0244797251582018D0 ]
    cf%ari = [ 0.0D0, -0.0321826755129339D0, 1.0D0, -0.0321826755129339D0, 0.0D0 ]
  ! one-sided weights, band-diagonal style
  ! no filter (current) or bad filter at boundary causes ringing in the solution
 !   cf%alb1(:,1,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
 !   cf%alb1(:,2,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%alb1(:,1,0) = [ 0.0D0 , 0.0D0  ,  1.375D0,   -0.75D0,    0.375D0 ]
    cf%alb1(:,2,0) = [ 0.0D0 , 0.3186803178523657D0, 0.2982740132694008D0, 0.3186803178523657D0, 0.0D0 ]
    cf%alb1(:,3,0) = cf%ali
    cf%alb2(:,1,0) = cf%alb1(ncl:1:-1,3,0)
    cf%alb2(:,2,0) = cf%alb1(ncl:1:-1,2,0)
    cf%alb2(:,3,0) = cf%alb1(ncl:1:-1,1,0)
    cf%arb1(:,1,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%arb1(:,2,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%arb1(:,2,0) = cf%ari
    cf%arb1(:,3,0) = cf%ari
    cf%arb2(:,1,0) = cf%arb1(ncr:1:-1,3,0)
    cf%arb2(:,2,0) = cf%arb1(ncr:1:-1,2,0)
    cf%arb2(:,3,0) = cf%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(cf%alb1(:,:,+1),cf%arb1(:,:,+1),cf%ali,cf%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(cf%alb1(:,:,-1),cf%arb1(:,:,-1),cf%ali,cf%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(cf%alb2(:,:,+1),cf%arb2(:,:,+1),cf%ali,cf%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(cf%alb2(:,:,-1),cf%arb2(:,:,-1),cf%ali,cf%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data : same as one-sided to maintain invertibility
!    cf%alb1(:,1,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%alb1(:,2,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%alb1(:,1,2) = [ 0.0D0 , 0.0D0  ,  1.375D0,   -0.75D0,    0.375D0 ]
    cf%alb1(:,2,2) = [ 0.0D0 , 0.3186803178523657D0, 0.2982740132694008D0, 0.3186803178523657D0, 0.0D0 ]
    cf%alb1(:,3,2) = cf%ali
    cf%alb2(:,1,2) = cf%alb1(ncl:1:-1,3,2)
    cf%alb2(:,2,2) = cf%alb1(ncl:1:-1,2,2)
    cf%alb2(:,3,2) = cf%alb1(ncl:1:-1,1,2)
    cf%arb1(:,1,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%arb1(:,2,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%arb1(:,2,2) = cf%ari
    cf%arb1(:,3,2) = cf%ari
    cf%arb2(:,1,2) = cf%arb1(ncr:1:-1,3,2)
    cf%arb2(:,2,2) = cf%arb1(ncr:1:-1,2,2)
    cf%arb2(:,3,2) = cf%arb1(ncr:1:-1,1,2)
  end subroutine cfamrcf

  subroutine cfamrfc(cf)  ! setup compact fine-to-coarse for 3x refinement
    implicit none
    TYPE(compact_weight), intent(out) :: cf
    INTEGER(c_int), PARAMETER :: nol = 2, nor = 2, nir = 2
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 3, nbc2 = 3  ! ?
    LOGICAL(c_bool),PARAMETER :: implicit_op = .true.
    REAL(c_double), PARAMETER :: two=2.0_c_double
    INTEGER :: i,j,k
    allocate( cf%ali(2*nol+1) )
    allocate( cf%ari(2*nor+1) )
    allocate( cf%alb1(2*nol+1,nbc1,-1:2), cf%alb2(2*nol+1,nbc2,-1:2) )
    allocate( cf%arb1(2*nor+1,nbc1,-1:2), cf%arb2(2*nor+1,nbc2,-1:2) )
    cf%description = compact_descriptor('famrfc31','explicit',0,0)
    if( implicit_op ) cf%description%complicity = 'implicit'
    call cf%description%info()
    cf%nol = nol  ! half-width of lhs stencil
    cf%nor = nor  ! half-width of outer rhs stencil
    cf%nir = nir  ! half-width of inner rhs stencil
    cf%ncl = ncl  ! total width of lhs stencil
    cf%ncr = ncr  ! total width of rhs stencil
    cf%nci = nci  ! parallel overlap of lhs stencil 
    cf%nst = nst  ! shift for staggered solution
    cf%nbc1 = nbc1  ! number of custom lower boundary points
    cf%nbc2 = nbc2  ! number of custom upper boundary points
    cf%null_option = 1 ! copy if null_op
    cf%implicit_op = implicit_op  ! lhs ignored if false
    if( nol == 0 ) cf%implicit_op = .false.
    cf%sh = 0.0d0  ! shift
  ! interior weights
    cf%ali = [ 0.0D0, -0.0321826755129339D0, 1.0D0, -0.0321826755129339D0, 0.0D0 ]
    cf%ari = [ 0.0244797251582018D0, 0.2207614172195584D0, 0.4451523642186118D0, 0.2207614172195584D0, 0.0244797251582018D0 ]
  ! one-sided weights, band-diagonal style
  ! no filter (current) or bad filter at boundary causes ringing in the solution
    cf%alb1(:,1,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%alb1(:,2,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%alb1(:,2,0) = cf%ali
    cf%alb1(:,3,0) = cf%ali
    cf%alb2(:,1,0) = cf%alb1(ncl:1:-1,3,0)
    cf%alb2(:,2,0) = cf%alb1(ncl:1:-1,2,0)
    cf%alb2(:,3,0) = cf%alb1(ncl:1:-1,1,0)
!    cf%arb1(:,1,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%arb1(:,2,0) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%arb1(:,1,0) = [ 0.0D0 , 0.0D0  ,  1.375D0,   -0.75D0,    0.375D0 ]
    cf%arb1(:,2,0) = [ 0.0D0 , 0.3186803178523657D0, 0.2982740132694008D0, 0.3186803178523657D0, 0.0D0 ]
    cf%arb1(:,3,0) = cf%ari
    cf%arb2(:,1,0) = cf%arb1(ncr:1:-1,3,0)
    cf%arb2(:,2,0) = cf%arb1(ncr:1:-1,2,0)
    cf%arb2(:,3,0) = cf%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(cf%alb1(:,:,+1),cf%arb1(:,:,+1),cf%ali,cf%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(cf%alb1(:,:,-1),cf%arb1(:,:,-1),cf%ali,cf%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(cf%alb2(:,:,+1),cf%arb2(:,:,+1),cf%ali,cf%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(cf%alb2(:,:,-1),cf%arb2(:,:,-1),cf%ali,cf%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data : same as one-sided to maintain invertibility
    cf%alb1(:,1,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%alb1(:,2,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%alb1(:,2,2) = cf%ali
    cf%alb1(:,3,2) = cf%ali
    cf%alb2(:,1,2) = cf%alb1(ncl:1:-1,3,2)
    cf%alb2(:,2,2) = cf%alb1(ncl:1:-1,2,2)
    cf%alb2(:,3,2) = cf%alb1(ncl:1:-1,1,2)
!    cf%arb1(:,1,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
!    cf%arb1(:,2,2) = [ 0.0D0 , 0.0D0  , 1.0D0 , 0.0D0  , 0.0D0 ]
    cf%arb1(:,1,2) = [ 0.0D0 , 0.0D0  ,  1.375D0,   -0.75D0,    0.375D0 ]
    cf%arb1(:,2,2) = [ 0.0D0 , 0.3186803178523657D0, 0.2982740132694008D0, 0.3186803178523657D0, 0.0D0 ]
    cf%arb1(:,3,2) = cf%ari
    cf%arb2(:,1,2) = cf%arb1(ncr:1:-1,3,2)
    cf%arb2(:,2,2) = cf%arb1(ncr:1:-1,2,2)
    cf%arb2(:,3,2) = cf%arb1(ncr:1:-1,1,2)
  end subroutine cfamrfc

  subroutine e4d4(ff)  ! setup 4th derivative, explicit
    implicit none
    TYPE(compact_weight), intent(out) :: ff
    INTEGER(c_int), PARAMETER :: nol = 0, nor = 3, nir = 3
    INTEGER(c_int), PARAMETER :: ncl=2*nol+1, ncr=nor+nir+1, nci=2*nol
    INTEGER(c_int), PARAMETER :: nst = 0, nbc1 = 3, nbc2 = 3
    LOGICAL(c_bool),PARAMETER :: implicit_op = .false.
    ! Fourth order 4th derivative
    DOUBLE PRECISION, PARAMETER :: agau = 28.0D0/3.0D0
    DOUBLE PRECISION, PARAMETER :: bgau = -6.5D0
    DOUBLE PRECISION, PARAMETER :: cgau = 2.0D0
    DOUBLE PRECISION, PARAMETER :: dgau = -1.0D0/6.0D0
    ! Second order 4th derivative
    DOUBLE PRECISION, PARAMETER :: a3gau = 6.0D0
    DOUBLE PRECISION, PARAMETER :: b3gau = -4.0D0
    DOUBLE PRECISION, PARAMETER :: c3gau = 1.0D0

    INTEGER :: i,j,k
    ff%description = compact_descriptor('d4exp','explicit',0,0)
    if( implicit_op ) ff%description%complicity = 'implicit'
    if( verbose) call ff%description%info()
    allocate( ff%ali(2*nol+1) )
    allocate( ff%ari(2*nor+1) )
    allocate( ff%alb1(2*nol+1,nbc1,-1:2), ff%alb2(2*nol+1,nbc2,-1:2) )
    allocate( ff%arb1(2*nor+1,nbc1,-1:2), ff%arb2(2*nor+1,nbc2,-1:2) )
    ff%nol = nol
    ff%nor = nor
    ff%nir = nir
    ff%ncl = ncl
    ff%ncr = ncr
    ff%nci = nci
    ff%nst = nst
    ff%nbc1 = nbc1
    ff%nbc2 = nbc2
    ff%null_option = 1 ! copy if null_op
    ff%implicit_op = implicit_op
    if( nol == 0 ) ff%implicit_op = .false.
    ff%sh = 0.0d0  ! shift
  ! interior weights
    !    ff%ali = 1.0D0
    !    ff%ari = [ dgau , cgau , bgau , agau , bgau , cgau , dgau ]
    ff%ali = 1.0D0
    ff%ari = [ dgau , cgau , bgau , agau , bgau , cgau , dgau ]
  ! one-sided weights, band-diagonal style
    ff%alb1(:,1,0) = 1.0D0
    ff%alb1(:,2,0) = 1.0D0
    ff%alb1(:,3,0) = 1.0D0
    ff%alb2(:,1,0) = ff%alb1(ncl:1:-1,3,0)
    ff%alb2(:,2,0) = ff%alb1(ncl:1:-1,2,0)
    ff%alb2(:,3,0) = ff%alb1(ncl:1:-1,1,0)
    ff%arb1(:,1,0) = [ 0.0D0   , 0.0D0    , 0.0D0    , 0.0D0    , 0.0D0  , 0.0D0     , 0.0D0  ]
    ff%arb1(:,2,0) = [ 0.0D0   , 0.0D0    , 0.0D0    , 0.0D0    , 0.0D0  , 0.0D0     , 0.0D0  ]
    ff%arb1(:,3,0) = [ 0.0D0   , c3gau    , b3gau    , a3gau    , b3gau  , c3gau     , 0.0D0  ]
    ff%arb2(:,1,0) = ff%arb1(ncr:1:-1,3,0)
    ff%arb2(:,2,0) = ff%arb1(ncr:1:-1,2,0)
    ff%arb2(:,3,0) = ff%arb1(ncr:1:-1,1,0)
  ! symmetric weights
    call lower_symm_weights(ff%alb1(:,:,+1),ff%arb1(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call lower_symm_weights(ff%alb1(:,:,-1),ff%arb1(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
    call upper_symm_weights(ff%alb2(:,:,+1),ff%arb2(:,:,+1),ff%ali,ff%ari,nol,nor,0,0,+1,+1)
    call upper_symm_weights(ff%alb2(:,:,-1),ff%arb2(:,:,-1),ff%ali,ff%ari,nol,nor,0,0,-1,-1)
  ! with extended boundary data
    ff%alb1(:,1,2) = 1.0D0
    ff%alb1(:,2,2) = 1.0D0
    ff%alb1(:,3,2) = 1.0D0
    ff%alb2(:,1,2) = ff%alb1(ncl:1:-1,3,2)
    ff%alb2(:,2,2) = ff%alb1(ncl:1:-1,2,2)
    ff%alb2(:,3,2) = ff%alb1(ncl:1:-1,1,2)
    ff%arb1(:,1,2) = ff%ari
    ff%arb1(:,2,2) = ff%ari
    ff%arb1(:,3,2) = ff%ari
    ff%arb2(:,1,2) = ff%arb1(ncr:1:-1,3,2)
    ff%arb2(:,2,2) = ff%arb1(ncr:1:-1,2,2)
    ff%arb2(:,3,2) = ff%arb1(ncr:1:-1,1,2)
  end subroutine e4d4

  

  subroutine lower_symm_weights_int(alb,arb,ali,ari,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:), intent(inout) :: alb,arb
    real(kind=c_double), dimension(:), intent(in) :: ali,ari
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,k,n,i1,i2
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
    do i=1,n
      alb(:,i) = ali
    end do
    do j=1,nol
      k = nol-j+1   ! nol:1:-1
      do i=1,k
        i1 = k-i+1  ! k:1:-1
        i2 = k+i+shl
        alb(i2,j) = alb(i2,j)+real(syml,kind=c_double)*alb(i1,j) ; alb(i1,j) = zero
      end do
    end do
    n = size(arb,2)
    do i=1,n
      arb(:,i) = ari
    end do
    do j=1,nor
      k = nor-j+1   ! nor:1:-1
      do i=1,k
        i1 = k-i+1  ! k:1:-1
        i2 = k+i+shr
        arb(i2,j) = arb(i2,j)+real(symr,kind=c_double)*arb(i1,j) ; arb(i1,j) = zero
      end do
    end do
  end subroutine lower_symm_weights_int
 
  subroutine upper_symm_weights_int(alb,arb,ali,ari,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:), intent(inout) :: alb,arb
    real(kind=c_double), dimension(:), intent(in) :: ali,ari
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,jj,k,n,noff
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
    do i=1,n
      alb(:,i) = ali
    end do
    do j=1,nol
      jj = n+1-j
      k = nol-j+1
      do i=1,k
        noff = size(ali)-k
        alb(noff-i+1-shl,jj) = alb(noff-i+1-shl,jj)+real(syml,kind=c_double)*alb(noff+i,jj) ; alb(noff+i,jj) = zero
      end do
    end do
    n = size(arb,2)
    do i=1,n
      arb(:,i) = ari
    end do
    do j=n-nor+1,n
      k = j-n+nor
      do i=1,k
        noff = size(ari)-k
        arb(noff-i+1-shr,j) = arb(noff-i+1-shr,j)+real(symr,kind=c_double)*arb(noff+i,j) ; arb(noff+i,j) = zero
      end do
    end do
  end subroutine upper_symm_weights_int

  subroutine lower_symm_weights_gen(alb,arb,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:), intent(inout) :: alb,arb
!    real(kind=c_double), dimension(:,:), intent(in) :: ali,ari
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,k,n
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
!    do i=1,n
!      alb(:,i) = ali(:,i)
!    end do
    do j=1,nol
      k = nol-j+1
      do i=1,k
        alb(2*k-i+1+shl,j) = alb(2*k-i+1+shl,j)+real(syml,kind=c_double)*alb(i,j) ; alb(i,j) = zero
      end do
    end do
    n = size(arb,2)
!    do i=1,n
!      arb(:,i) = ari(:,i)
!    end do
    do j=1,nor
      k = nor-j+1
      do i=1,k
        arb(2*k-i+1+shr,j) = arb(2*k-i+1+shr,j)+real(symr,kind=c_double)*arb(i,j) ; arb(i,j) = zero
      end do
    end do
  end subroutine lower_symm_weights_gen
 
  subroutine upper_symm_weights_gen(alb,arb,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:), intent(inout) :: alb,arb
!    real(kind=c_double), dimension(:,:), intent(in) :: ali,ari
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,k,n,noff
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
!    do i=1,n
!      alb(:,i) = ali(:,i)
!    end do
    do j=n-nol+1,n
      k = j-n+nol
      do i=1,k
        noff = size(alb,1)-k
        alb(noff-i+1-shl,j) = alb(noff-i+1-shl,j)+real(syml,kind=c_double)*alb(noff+i,j) ; alb(noff+i,j) = zero
      end do
    end do
    n = size(arb,2)
!    do i=1,n
!      arb(:,i) = ari(:,i)
!    end do
    do j=n-nor+1,n
      k = j-n+nor
      do i=1,k
        noff = size(arb,1)-k
        arb(noff-i+1-shr,j) = arb(noff-i+1-shr,j)+real(symr,kind=c_double)*arb(noff+i,j) ; arb(noff+i,j) = zero
      end do
    end do
  end subroutine upper_symm_weights_gen

  subroutine lower_symm_weights_both(alb,arb,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:,-1:), intent(inout) :: alb,arb
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,k,n
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
    do j=1,nol
      k = nol-j+1
      do i=1,k
        alb(2*k-i+1+shl,j,-1) = alb(2*k-i+1+shl,j,-1)-real(syml,kind=c_double)*alb(i,j,-1) ; alb(i,j,-1) = zero
        alb(2*k-i+1+shl,j,+1) = alb(2*k-i+1+shl,j,+1)+real(syml,kind=c_double)*alb(i,j,+1) ; alb(i,j,+1) = zero
      end do
    end do
    n = size(arb,2)
    do j=1,nor
      k = nor-j+1
      do i=1,k
        arb(2*k-i+1+shr,j,-1) = arb(2*k-i+1+shr,j,-1)-real(symr,kind=c_double)*arb(i,j,-1) ; arb(i,j,-1) = zero
        arb(2*k-i+1+shr,j,+1) = arb(2*k-i+1+shr,j,+1)+real(symr,kind=c_double)*arb(i,j,+1) ; arb(i,j,+1) = zero
      end do
    end do
  end subroutine lower_symm_weights_both
 
  subroutine upper_symm_weights_both(alb,arb,nol,nor,shl,shr,syml,symr)
    implicit none
    real(kind=c_double), dimension(:,:,-1:), intent(inout) :: alb,arb
    INTEGER(c_int), intent(in) :: nol,nor,shl,shr,symr,syml
    INTEGER :: i,j,k,n,noff
    REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
    n = size(alb,2)
    do j=n-nol+1,n
      k = j-n+nol
      do i=1,k
        noff = size(alb,1)-k
        alb(noff-i+1-shl,j,-1) = alb(noff-i+1-shl,j,-1)-real(syml,kind=c_double)*alb(noff+i,j,-1) ; alb(noff+i,j,-1) = zero
        alb(noff-i+1-shl,j,+1) = alb(noff-i+1-shl,j,+1)+real(syml,kind=c_double)*alb(noff+i,j,+1) ; alb(noff+i,j,+1) = zero
      end do
    end do
    n = size(arb,2)
    do j=n-nor+1,n
      k = j-n+nor
      do i=1,k
        noff = size(arb,1)-k
        arb(noff-i+1-shr,j,-1) = arb(noff-i+1-shr,j,-1)-real(symr,kind=c_double)*arb(noff+i,j,-1) ; arb(noff+i,j,-1) = zero
        arb(noff-i+1-shr,j,+1) = arb(noff-i+1-shr,j,+1)+real(symr,kind=c_double)*arb(noff+i,j,+1) ; arb(noff+i,j,+1) = zero
      end do
    end do
  end subroutine upper_symm_weights_both

end module LES_stencils

