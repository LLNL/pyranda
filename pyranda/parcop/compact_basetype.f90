module LES_compact_basetype
  USE iso_c_binding
  USE MPI
  !USE LES_input, ONLY : directcom
  USE LES_stencils
  USE LES_comm,  ONLY : mpierr
  use LES_pentadiagonal, ONLY : bpentLUD1, ppentLUD1, bpentLUS1, ppentLUS1, &
    btrid_block4_lud, ptrid_block4_lud

  IMPLICIT NONE

  integer :: directcom = 1
  
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double, one=1.0_c_double
  LOGICAL(c_bool) :: debug=.false.

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
    real(kind=c_double), dimension(:,:), pointer :: art => null() ! transposed rhs (matrix compat.)
    real(kind=c_double), dimension(:,:), pointer :: ar  => null() ! rhs
    real(kind=c_double), dimension(:,:), pointer :: al  => null() ! lhs
    real(kind=c_double), dimension(:,:), pointer :: sm  => null() ! lhs^1*rhs
    real(kind=c_double), dimension(:,:), pointer :: smt => null() ! lhs^1*rhs transpose
    real(kind=c_double), dimension(:,:), pointer :: rc  => null() ! rhs overlap 
    real(kind=c_double), dimension(:,:,:,:), pointer :: aa => null() ! parallel weights
  contains   ! generic
    procedure :: setup => setup_compact_op1
    procedure :: remove => remove_compact_op1
  END TYPE compact_op1

  TYPE mesh1_type  ! 1D
    INTEGER(c_int) :: m,n                    ! mesh size  -> a[x,y,z], n[x,y,z]
    real(kind=c_double) :: d,bv1,bvn         ! grid spacing/metric -> d[x,y,z], [x,y,z][1,n]
    CHARACTER(KIND=c_char,LEN=4) :: bc1,bcn  ! boundary conditions (string) -> b[x,y,z][1,n]
  END TYPE mesh1_type
  
  TYPE comm1_type  ! 1D ! for testing
    logical(c_bool) :: periodic  ! -> periodic[x,y,z]
    INTEGER(c_int) :: hash,np,id,lo,hi,range(2)  ! -> [x,y,z]com,com_{np,id,lo,hi},range
  END TYPE comm1_type
  
! parallel overlap (mpi_allgather)
  real(kind=8), dimension(:,:,:), allocatable :: dvopx,dvopy,dvopz
  real(kind=8), dimension(:,:,:,:), allocatable :: dvox,dvoy,dvoz
! halo buffer (mpi_sendrecv)
  real(kind=8), dimension(:,:,:), allocatable :: vbr1x,vbr2x,vbs1x,vbs2x
  real(kind=8), dimension(:,:,:), allocatable :: vbr1y,vbr2y,vbs1y,vbs2y
  real(kind=8), dimension(:,:,:), allocatable :: vbr1z,vbr2z,vbs1z,vbs2z

contains

  subroutine setup_compact_op1(op,wgt,com,msh,bc,null_op)  ! direction agnostic
    IMPLICIT NONE
    class(compact_op1) :: op
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
    ! free pre-existing setup
    call op%remove()
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
    allocate( op%sm(m,1:m+2*op%nor) ) ; op%sm = zero
    allocate( op%smt(1:m+2*op%nor,m) ) ; op%smt = zero
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
    do j=1,m ; op%sm(j,j:j+2*op%nor) = op%art(j,:) ; end do  ! rhs matrix

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
      call mpi_allgather(rop,ni*ni,MPI_DOUBLE_PRECISION,ro,ni*ni,MPI_DOUBLE_PRECISION,op%hash,mpierr)
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
    if( np==1 .and. op%periodic ) then
      do j=1,m+2*op%nor ; call ppentLUS1(op%al, op%sm(:,j), m) ; end do  ! lhs^1*rhs matrix
    else
      do j=1,m+2*op%nor ; call bpentLUS1(op%al, op%sm(:,j), m) ; end do  ! lhs^1*rhs matrix
    endif
    endif ! implicit
    op%smt = transpose(op%sm)
    deallocate( gal, gar )
    deallocate( ro )
    deallocate( rop, al1, al2 )
  end subroutine setup_compact_op1

  subroutine remove_compact_op1(op)  ! free pointer arrays manually
    class(compact_op1) :: op
    if( associated(op%art) ) deallocate(op%art) ; op%art => null()
    if( associated(op%ar ) ) deallocate(op%ar)  ; op%ar  => null()
    if( associated(op%al ) ) deallocate(op%al)  ; op%al  => null()
    if( associated(op%sm ) ) deallocate(op%sm)  ; op%sm  => null()
    if( associated(op%smt) ) deallocate(op%smt) ; op%smt => null()
    if( associated(op%rc ) ) deallocate(op%rc)  ; op%rc  => null()
    if( associated(op%aa ) ) deallocate(op%aa)  ; op%aa  => null()
  end subroutine remove_compact_op1

end module LES_compact_basetype

