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
!===================================================================================================
 MODULE LES_comm
!===================================================================================================
  USE iso_c_binding
  !USE MPI
  !USE LES_input, ONLY : xop,yop,zop,xdice,ydice,zdice,compressible,calcSpectrum,print_screen
  USE LES_patch, ONLY : patch_type
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER(c_int), PARAMETER :: master = 0         ! master process
  INTEGER(c_int)            :: LES_comm_world = 0 ! settable world comm       *APIVAR*
  INTEGER(c_int)            :: world_id = 0       ! world-rank of MPI process *APIVAR*
  INTEGER(c_int)            :: world_np = 1       ! number of MPI processes
  INTEGER(c_int), DIMENSION(MPI_STATUS_SIZE) :: mpistatus
  INTEGER(c_int) :: mpisource,mpitag,mpierr

  LOGICAL(c_bool), PARAMETER :: compressible = .TRUE.
  LOGICAL(c_bool), PARAMETER :: calcSpectrum = .FALSE.
  LOGICAL(c_bool), PARAMETER :: print_screen = .TRUE.
  INTEGER(c_int) :: xop = 2                   ! index upon which x derivatives and filters are applied
  INTEGER(c_int) :: yop = 2                   ! index upon which y derivatives and filters are applied
  INTEGER(c_int) :: zop = 3                   ! index upon which z derivatives and filters are applied
  INTEGER(c_int) :: xdice = 3                 ! index/plane subdivided in x-transpose
  INTEGER(c_int) :: ydice = 3                 ! index/plane subdivided in y-transpose
  INTEGER(c_int) :: zdice = 1                 ! index/plane subdivided in z-transpose
  
  TYPE comm_type
    ! PATCH COMMUNICATOR ---------------------------------------------------------------------------
    INTEGER(c_int) :: patcom,patcom_np,patcom_id
    ! CARTESIAN COMMUNICATORS ----------------------------------------------------------------------
    INTEGER(c_int) :: xyzcom,xyzcom_np,xyzcom_id
    INTEGER(c_int) :: xycom,xycom_np,xycom_id
    INTEGER(c_int) :: xzcom,xzcom_np,xzcom_id
    INTEGER(c_int) :: yzcom,yzcom_np,yzcom_id
    INTEGER(c_int) :: xcom,xcom_np,xcom_id,xcom_lo,xcom_hi
    INTEGER(c_int) :: ycom,ycom_np,ycom_id,ycom_lo,ycom_hi
    INTEGER(c_int) :: zcom,zcom_np,zcom_id,zcom_lo,zcom_hi
    INTEGER(c_int), DIMENSION(2) :: xrange,yrange,zrange
    ! TRANSPOSES -----------------------------------------------------------------------------------
    INTEGER(c_int) :: bx_xtran ! number of x grid points per processor after x transpose
    INTEGER(c_int) :: by_xtran ! number of y grid points per processor after x transpose
    INTEGER(c_int) :: bz_xtran ! number of z grid points per processor after x transpose
    INTEGER(c_int) :: bx_ytran ! number of x grid points per processor after y transpose
    INTEGER(c_int) :: by_ytran ! number of y grid points per processor after y transpose
    INTEGER(c_int) :: bz_ytran ! number of z grid points per processor after y transpose
    INTEGER(c_int) :: bx_ztran ! number of x grid points per processor after z transpose
    INTEGER(c_int) :: by_ztran ! number of y grid points per processor after z transpose
    INTEGER(c_int) :: bz_ztran ! number of z grid points per processor after z transpose
    CONTAINS
     PROCEDURE :: setup => setup_comm
     PROCEDURE :: remove => remove_comm
!     FINAL :: remove_mesh     	
  END TYPE comm_type
  
!===================================================================================================
  CONTAINS
!=================================================================================================== 

   SUBROUTINE setup_comm(comm_data,patch_data)
    IMPLICIT NONE
    CLASS(comm_type),  INTENT(OUT) :: comm_data ! Auto deallocation of type components
    CLASS(patch_type), INTENT(IN)  :: patch_data
    LOGICAL, PARAMETER :: reorder = .TRUE.
    INTEGER, PARAMETER :: ndim = 3
    CHARACTER(KIND=c_char,LEN=3), PARAMETER :: fileorder = 'XYZ' ! processor-ordering restart files (first coordinate varies the fastest) 
    INTEGER, DIMENSION(ndim), PARAMETER :: pdir = (/ 3, 2, 1 /)  ! goes with fileorder
    INTEGER :: xyzmap(3),xymap(2),xzmap(2),yzmap(2)              ! maps
    INTEGER :: xdice_valid(3),bxdim(3),bxdim_valid(3,3)          ! register of valid x-transposes
    INTEGER :: ydice_valid(3),bydim(3),bydim_valid(3,3)          ! register of valid y-transposes
    INTEGER :: zdice_valid(3),bzdim(3),bzdim_valid(3,3)          ! register of valid z-transposes
    INTEGER, DIMENSION(ndim) :: pdim,pmap,bdim
    LOGICAL, DIMENSION(ndim) :: periodic,remain
    INTEGER :: xdirection,ydirection,zdirection,i,procnum,cf,pf,fatal_err
    LOGICAL :: invalid,unset
     ASSOCIATE(color     => patch_data%color,     key       => patch_data%key,                                          &
               nx        => patch_data%nx,        ny        => patch_data%ny,        nz        => patch_data%nz,        &
               px        => patch_data%px,        py        => patch_data%py,        pz        => patch_data%pz,        &
               ax        => patch_data%ax,        ay        => patch_data%ay,        az        => patch_data%az,        &
               bx1       => patch_data%bx1,       by1       => patch_data%by1,       bz1       => patch_data%bz1,       &
               periodicx => patch_data%periodicx, periodicy => patch_data%periodicy, periodicz => patch_data%periodicz, &
               xrange    => comm_data%xrange,     yrange    => comm_data%yrange,     zrange    => comm_data%zrange,     &
               patcom    => comm_data%patcom,     patcom_np => comm_data%patcom_np,  patcom_id => comm_data%patcom_id,  &
               xyzcom    => comm_data%xyzcom,     xyzcom_np => comm_data%xyzcom_np,  xyzcom_id => comm_data%xyzcom_id,  &
               xycom     => comm_data%xycom,      xycom_np  => comm_data%xycom_np,   xycom_id  => comm_data%xycom_id,   &
               xzcom     => comm_data%xzcom,      xzcom_np  => comm_data%xzcom_np,   xzcom_id  => comm_data%xzcom_id,   &
               yzcom     => comm_data%yzcom,      yzcom_np  => comm_data%yzcom_np,   yzcom_id  => comm_data%yzcom_id,   &
               xcom      => comm_data%xcom,       xcom_np   => comm_data%xcom_np,    xcom_id   => comm_data%xcom_id,    &
	       xcom_lo   => comm_data%xcom_lo,    xcom_hi   => comm_data%xcom_hi,                                       &
               ycom      => comm_data%ycom,       ycom_np   => comm_data%ycom_np,    ycom_id   => comm_data%ycom_id,    &
	       ycom_lo   => comm_data%ycom_lo,    ycom_hi   => comm_data%ycom_hi,                                       &
               zcom      => comm_data%zcom,       zcom_np   => comm_data%zcom_np,    zcom_id   => comm_data%zcom_id,    &
	       zcom_lo   => comm_data%zcom_lo,    zcom_hi   => comm_data%zcom_hi,                                       &
               bx_xtran  => comm_data%bx_xtran,   by_xtran  => comm_data%by_xtran,   bz_xtran  => comm_data%bz_xtran,   &
               bx_ytran  => comm_data%bx_ytran,   by_ytran  => comm_data%by_ytran,   bz_ytran  => comm_data%bz_ytran,   &
               bx_ztran  => comm_data%bx_ztran,   by_ztran  => comm_data%by_ztran,   bz_ztran  => comm_data%bz_ztran)
! patch communicator
     CALL MPI_COMM_SPLIT(LES_comm_world, color, key, patcom, mpierr);
     CALL MPI_COMM_SIZE(patcom,patcom_np,mpierr)	
     CALL MPI_COMM_RANK(patcom,patcom_id,mpierr)       
! user-specified ordering of processor layout
!     SELECT CASE(fileorder) ! highest number varies fastest
!     CASE('ZYX')
!      pdir = (/ 1, 2, 3 /)
!     CASE('XYZ')
!      pdir = (/ 3, 2, 1 /)
!     CASE('ZXY')
!      pdir = (/ 2, 1, 3 /)
!     CASE DEFAULT
!      IF (patcom_id == master) PRINT *, 'fileorder:',fileorder,' not supported'
!      CALL MPI_FINALIZE(mpierr)
!      STOP
!     END SELECT
!     DO procnum=1,ndim
!      IF (.NOT. ANY(pdir == procnum)) THEN
!       IF (patcom_id == master) PRINT *, 'Bad process ordering. Using default.'
!       pdir = (/ (i,i=1,ndim) /)
!       EXIT
!      END IF
!     END DO
! user-specified processor dimensions
     pdim(pdir) = (/ px, py, pz /)
     IF (PRODUCT(pdim) > patcom_np) THEN
      IF (patcom_id == master) PRINT *, 'Too few processors available: ',patcom_np
      CALL MPI_FINALIZE(mpierr)
      STOP
     END IF
     IF ((PRODUCT(pdim) < patcom_np) .AND. (patcom_id == master)) THEN
      PRINT *, 'Processor underutilization detected: patcom_np = ',patcom_np,' pdim = ',pdim
      CALL MPI_FINALIZE(mpierr)
      STOP
     END IF
!     ax = nx/pdim(pdir(1))
!     ay = ny/pdim(pdir(2))
!     az = nz/pdim(pdir(3))
     periodic(pdir) = (/ periodicx, periodicy, periodicz /)
! xyzcom
     CALL MPI_CART_CREATE(patcom,ndim,pdim,periodic,reorder,xyzcom,mpierr)
     CALL MPI_COMM_RANK(xyzcom,xyzcom_id,mpierr)
     CALL MPI_COMM_SIZE(xyzcom,xyzcom_np,mpierr)
     CALL MPI_CART_COORDS(xyzcom,xyzcom_id,ndim,xyzmap,mpierr)
! xycom
     remain(pdir) = (/ .TRUE., .TRUE., .FALSE. /)
     CALL MPI_CART_SUB(xyzcom,remain,xycom,mpierr)
     CALL MPI_COMM_RANK(xycom,xycom_id,mpierr)
     CALL MPI_COMM_SIZE(xycom,xycom_np,mpierr)
     CALL MPI_CART_COORDS(xycom,xycom_id,2,xymap,mpierr)
! xzcom
     remain(pdir) = (/ .TRUE., .FALSE., .TRUE. /)
     CALL MPI_CART_SUB(xyzcom,remain,xzcom,mpierr)
     CALL MPI_COMM_RANK(xzcom,xzcom_id,mpierr)
     CALL MPI_COMM_SIZE(xzcom,xzcom_np,mpierr)
     CALL MPI_CART_COORDS(xzcom,xzcom_id,2,xzmap,mpierr)
! yzcom
     remain(pdir) = (/ .FALSE., .TRUE., .TRUE. /)
     CALL MPI_CART_SUB(xyzcom,remain,yzcom,mpierr)
     CALL MPI_COMM_RANK(yzcom,yzcom_id,mpierr)
     CALL MPI_COMM_SIZE(yzcom,yzcom_np,mpierr)
     CALL MPI_CART_COORDS(yzcom,yzcom_id,2,yzmap,mpierr)
! xcom
     remain(pdir) = (/ .TRUE., .FALSE., .FALSE. /)
     CALL MPI_CART_SUB(xyzcom,remain,xcom,mpierr)
     CALL MPI_COMM_RANK(xcom,xcom_id,mpierr)
     CALL MPI_COMM_SIZE(xcom,xcom_np,mpierr)
     CALL MPI_CART_SHIFT(xcom,0,1,xcom_lo,xcom_hi,mpierr)
! ycom
     remain(pdir) = (/ .FALSE., .TRUE., .FALSE. /)
     CALL MPI_CART_SUB(xyzcom,remain,ycom,mpierr)
     CALL MPI_COMM_RANK(ycom,ycom_id,mpierr)
     CALL MPI_COMM_SIZE(ycom,ycom_np,mpierr)
     CALL MPI_CART_SHIFT(ycom,0,1,ycom_lo,ycom_hi,mpierr)
! zcom
     remain(pdir) = (/ .FALSE., .FALSE., .TRUE. /)
     CALL MPI_CART_SUB(xyzcom,remain,zcom,mpierr)
     CALL MPI_COMM_RANK(zcom,zcom_id,mpierr)
     CALL MPI_COMM_SIZE(zcom,zcom_np,mpierr)
     CALL MPI_CART_SHIFT(zcom,0,1,zcom_lo,zcom_hi,mpierr)
! Sanity checks
     if (xyzcom_id == master) then
       ! These will be the same on all procs, just check on master
       fatal_err = 0
       if (pdim(pdir(1))*ax .ne. nx) then
          print *, 'Bad size in X direction: nx,px,ax = ',nx,pdim(pdir(1)),ax
          fatal_err = 1
       end if
       if (pdim(pdir(2))*ay .ne. ny) then
          print *, 'Bad size in Y direction: ny,py,ay = ',ny,pdim(pdir(2)),ay
          fatal_err = 1
       end if
       if (pdim(pdir(3))*az .ne. nz) then
          print *, 'Bad size in Z direction: nz,pz,az = ',nz,pdim(pdir(3)),az
          fatal_err = 1
       end if
       if (fatal_err == 1) then
          print *, 'Bad input, aborting run...'
          call MPI_ABORT(xyzcom,1,mpierr)
       end if
     end if
! coordinates
     xrange(1) =  xyzmap(pdir(1))*ax+1
     xrange(2) = (xyzmap(pdir(1))+1)*ax
     yrange(1) =  xyzmap(pdir(2))*ay+1
     yrange(2) = (xyzmap(pdir(2))+1)*ay
     zrange(1) =  xyzmap(pdir(3))*az+1
     zrange(2) = (xyzmap(pdir(3))+1)*az
! ensure nondistributed dimensions are not transposed
     if( pdim(pdir(1)) == 1 ) xop = 1
     if( pdim(pdir(2)) == 1 ) yop = 2
     if( pdim(pdir(3)) == 1 ) zop = 3
! complex factor for Poisson solve
     cf = 1
     if( .not. compressible .and. any( periodic ) .and. product(pdim) > 1 ) cf = 2
! allowed x transpose decompositions -- this is the default order
     pf = cf*px
     if( px == 1 ) pf=px
     xdice_valid = 0
     bxdim_valid = 0
     if( mod(az,pf) == 0 ) then
      bxdim_valid(1,3) = ax*px
      bxdim_valid(2,3) = ay
      bxdim_valid(3,3) = az/px
      xdice_valid(1) = 3
     endif
     if( mod(ay,pf) == 0 ) then
      bxdim_valid(1,2) = ax*px
      bxdim_valid(2,2) = ay/px
      bxdim_valid(3,2) = az
      xdice_valid(2) = 2
     endif
     if( mod(ay*az,pf) == 0 ) then
      bxdim_valid(1,1) = ax*px
      bxdim_valid(2,1) = ay*az/px
      bxdim_valid(3,1) = 1
      xdice_valid(3) = 1
     endif
!     IF( xyzcom_id == master ) then
!      print *,'Valid x transposes: xdice, dimensions for xop=1'
!      do i=1,3
!        if( xdice_valid(i) > 0 ) print *,xdice_valid(i),bxdim_valid(:,xdice_valid(i))
!      end do
!     endif
     if( all( xdice_valid == 0 ) .AND. (calcSpectrum .OR. (.NOT. compressible)) ) then
      IF ( xyzcom_id == master ) PRINT *, 'Warning: No x-transposes can be formed.'
      STOP
     end if
! allowed y transpose decompositions -- this is the default order
     pf = cf*py
     if( py == 1 ) pf=py
     ydice_valid = 0
     bydim_valid = 0
     if( mod(az,pf) == 0 ) then
      bydim_valid(2,3) = ay*py
      bydim_valid(1,3) = ax
      bydim_valid(3,3) = az/py
      ydice_valid(1) = 3
     endif
     if( mod(ax,pf) == 0 ) then
      bydim_valid(2,1) = ay*py
      bydim_valid(1,1) = ax/py
      bydim_valid(3,1) = az
      ydice_valid(2) = 1
     endif
     if( mod(ax*az,pf) == 0 ) then
      bydim_valid(2,2) = ay*py
      bydim_valid(3,2) = az*ax/py
      bydim_valid(1,2) = 1
      ydice_valid(3) = 2
     endif
!     IF( xyzcom_id == master ) then
!      print *,'Valid y transposes: ydice, dimensions for yop=2'
!      do i=1,3
!        if( ydice_valid(i) > 0 ) print *,ydice_valid(i),bydim_valid(:,ydice_valid(i))
!      end do
!     end if
     if( all( ydice_valid == 0 ) .AND. (calcSpectrum .OR. (.NOT. compressible)) ) then
      IF ( xyzcom_id == master ) PRINT *, 'Warning: No y-transposes can be formed.'
      STOP
     end if
! allowed z transpose decompositions -- this is the default order
     pf = cf*pz
     if( pz == 1 ) pf=pz
     zdice_valid = 0
     bzdim_valid = 0
     if( mod(ax,pf) == 0 ) then
      bzdim_valid(3,1) = az*pz
      bzdim_valid(1,1) = ax/pz
      bzdim_valid(2,1) = ay
      zdice_valid(1) = 1
     endif
     if( mod(ay,pf) == 0 ) then
      bzdim_valid(3,2) = az*pz
      bzdim_valid(1,2) = ax
      bzdim_valid(2,2) = ay/pz
      zdice_valid(2) = 2
     endif
     if( mod(ax*ay,pf) == 0 ) then
      bzdim_valid(3,3) = az*pz
      bzdim_valid(1,3) = ax*ay/pz
      bzdim_valid(2,3) = 1
      zdice_valid(3) = 3
     endif
!     IF( xyzcom_id == master ) then
!      print *,'Valid z transposes: zdice, dimensions for zop=3'
!      do i=1,3
!        if( zdice_valid(i) > 0 ) print *,zdice_valid(i),bzdim_valid(:,zdice_valid(i))
!      end do
!     end if
     if( all( zdice_valid == 0 ) .AND. (calcSpectrum .OR. (.NOT. compressible)) ) then
      IF ( xyzcom_id == master ) PRINT *, 'Warning: No z-transposes can be formed.'
      STOP
     end if
! check *dice prefs
     unset = ( xdice == 0 )
     if( unset ) then
      invalid = .false.
     else
      invalid = all( xdice_valid /= xdice )
     endif
     if( invalid .or. unset ) then  ! change xdice to first valid option
      do i=1,3
        if( xdice_valid(i) > 0 ) then
          xdice = xdice_valid(i)
          exit
        endif
      end do
      IF ( xyzcom_id == master ) then
!       if( invalid ) PRINT *, 'Warning: Requested xdice is invalid. Using ',xdice
!       if( unset ) print *, 'Using default xdice = ',xdice
      end if
     end if
     bxdim = cshift(bxdim_valid(:,xdice),1-xop)
     bx_xtran = bxdim(1)
     by_xtran = bxdim(2)
     bz_xtran = bxdim(3)
     unset = ( ydice == 0 )
     if( unset ) then
      invalid = .false.
     else
      invalid = all( ydice_valid /= ydice )
     endif
     if( invalid .or. unset ) then  ! change ydice to first valid option
      do i=1,3
        if( ydice_valid(i) > 0 ) then
          ydice = ydice_valid(i)
          exit
        endif
      end do
      IF ( xyzcom_id == master ) then
!       if( invalid ) PRINT *, 'Warning: Requested ydice is invalid. Using ',ydice
!       if( unset ) print *, 'Using default ydice = ',ydice
      end if
     endif
     bydim = cshift(bydim_valid(:,ydice),2-yop)
     bx_ytran = bydim(1)
     by_ytran = bydim(2)
     bz_ytran = bydim(3)
     unset = ( zdice == 0 )
     if( unset ) then
      invalid = .false.
     else
      invalid = all( zdice_valid /= zdice )
     endif
     if( invalid .or. unset ) then  ! change zdice to first valid option
      do i=1,3
        if( zdice_valid(i) > 0 ) then
          zdice = zdice_valid(i)
          exit
        endif
      end do
      IF ( xyzcom_id == master ) then
!       if( invalid ) PRINT *, 'Warning: Requested zdice is invalid. Using ',zdice
!       if( unset ) print *, 'Using default zdice = ',zdice
      end if
     end if
     bzdim = cshift(bzdim_valid(:,zdice),3-zop)
     bx_ztran = bzdim(1)
     by_ztran = bzdim(2)
     bz_ztran = bzdim(3)
! print processor maps etc.
!    IF (print_screen) THEN
!     SELECT CASE(xyzcom_id)
!     CASE(master)
!      WRITE(6,'(A,3I7)') ' GRID                 : ',nx,ny,nz
!      WRITE(6,'(A,3I7)') ' PROCESSOR ARRANGEMENT: ',px,py,pz
!      WRITE(6,'(A,3I7)') ' GRID PER PROCESSOR   : ',ax,ay,az
!     END SELECT
!    END IF
!    IF ( ((spectral) .OR. (.NOT. compressible)) .AND. (.NOT. ares_library) ) THEN
!    IF (writeVis .OR. writeGrid) THEN
!     WRITE (pmapFile,'(A)') './procmap'
!     IF ( xyzcom_id == master ) THEN
!      OPEN(10,FILE=TRIM(pmapFile),FORM='formatted',ACTION='write')
!      WRITE(10,'(A,3I7)') ' X-TRANSPOSED GRID/CPU: ',bx_xtran,by_xtran,bz_xtran
!      WRITE(10,'(A,3I7)') ' Y-TRANSPOSED GRID/CPU: ',bx_ytran,by_ytran,bz_ytran
!      WRITE(10,'(A,3I7)') ' Z-TRANSPOSED GRID/CPU: ',bx_ztran,by_ztran,bz_ztran
!      WRITE(10,'(A,3I2,A,3I2,A)') ' DEFAULT TRANSPOSE OPTIONS: iop=(',xop,yop,zop,' ), idice=(',xdice,ydice,zdice,' )'
!      WRITE(10,*) pdim
!      WRITE(10,'(A8,2x,3(A5,1x))') 'XYZ_RANK','X','Y','Z'
!      DO procnum=0,xyzcom_np-1
!        CALL MPI_CART_COORDS(xyzcom,procnum,ndim,pmap,mpierr)
!        WRITE(10,'(i8,2x,3(i5,1x))') procnum,pmap
!      END DO
!      CLOSE(10)
!     END IF
!    END IF
     END ASSOCIATE
   END SUBROUTINE setup_comm

   SUBROUTINE remove_comm(comm_data)
    CLASS(comm_type), INTENT(INOUT) :: comm_data
      IF (comm_data%zcom   .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%zcom,   mpierr)
      IF (comm_data%ycom   .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%ycom,   mpierr)
      IF (comm_data%xcom   .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%xcom,   mpierr)
      IF (comm_data%yzcom  .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%yzcom,  mpierr)
      IF (comm_data%xzcom  .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%xzcom,  mpierr)
      IF (comm_data%xycom  .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%xycom,  mpierr)
      IF (comm_data%xyzcom .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%xyzcom, mpierr)
      IF (comm_data%patcom .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE( comm_data%patcom, mpierr)
   END SUBROUTINE remove_comm
      
!===================================================================================================
 END MODULE LES_comm
!===================================================================================================
