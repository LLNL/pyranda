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
 MODULE LES_patch ! IJK GRID SPECIFICATIONS
!===================================================================================================  
  USE iso_c_binding
  IMPLICIT NONE
  
  TYPE patch_type                                                   ! color and key are for MPI_COMM_SPLIT
    INTEGER(c_int)                 :: color                         ! Processes with same color will be placed in same communicator.
    INTEGER(c_int)                 :: key                           ! For same colors, same keys preserve ranks, whereas different keys indicate rank order
    INTEGER(c_int)                 :: coordsys                      ! coordinate system: 0=Cartesian, 1=Cylindrical, 2=Spherical, 3=Curvilinear
    INTEGER(c_int)                 :: nx,ny,nz                      ! dimensions 
    INTEGER(c_int)                 :: px,py,pz                      ! processor grid
    REAL(c_double)                 :: x1,xn,y1,yn,z1,zn             ! bounding box (cell faces or nodes)
    CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn       ! boundary conditions
    REAL(c_double)                 :: simtime                       ! simulation time
    REAL(c_double)                 :: dt = 0.0D0                    ! n   time increment
    REAL(c_double)                 :: dtold   = -1.0D0              ! n-1 time increment
    REAL(c_double)                 :: dtolder = -1.0D0              ! n-2 time increment
    CHARACTER(KIND=c_char, LEN=10) :: stability = '  startup '      ! stability indicator
    INTEGER(c_int), DIMENSION(3)   :: ijkstab = [0,0,0]             ! indices of stability determining node
    INTEGER(c_int)                 :: ax,ay,az                      ! per-processor dimensions
    REAL(c_double)                 :: dx,dy,dz                      ! distances between cell centers for Cartesian meshes
    INTEGER(c_int)                 :: isymX = 1                     ! multiplier for symmetric and anti-symmetric x boundaries
    INTEGER(c_int)                 :: isymY = 1                     ! multiplier for symmetric and anti-symmetric y boundaries
    INTEGER(c_int)                 :: isymZ = 1                     ! multiplier for symmetric and anti-symmetric z boundaries
    LOGICAL                        :: periodicx,periodicy,periodicz ! periodic directions (FORTRAN LOGICALS ARE 4 BYTE OBJECTS!)
!    LOGICAL(c_bool)                :: periodicx,periodicy,periodicz ! periodic directions (C LOGICALS ARE 1 BYTE OBJECTS! THIS FOOLS MPI, WHICH AUTOMATICALLY CONVERTS FROM FORTRAN TO C.))
    CONTAINS
     PROCEDURE :: setup => setup_patch
     PROCEDURE :: remove => remove_patch
!     FINAL :: remove_patch    	
  END TYPE patch_type
  
!===================================================================================================
  CONTAINS
!=================================================================================================== 

   SUBROUTINE setup_patch(patch_data,color,key,coordsys,nx,ny,nz,px,py,pz,x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn,simtime)
    IMPLICIT NONE
    CLASS(patch_type),            INTENT(OUT) :: patch_data ! Auto deallocation of all type components
    INTEGER(c_int),               INTENT(IN)  :: color,key,coordsys,nx,ny,nz,px,py,pz
    REAL(c_double),               INTENT(IN)  :: x1,xn,y1,yn,z1,zn
    CHARACTER(KIND=c_char,LEN=*), INTENT(IN)  :: bx1,bxn,by1,byn,bz1,bzn
    REAL(c_double),               INTENT(IN)  :: simtime
     patch_data%color = color
     patch_data%key = key
     patch_data%coordsys = coordsys
     patch_data%nx = nx
     patch_data%ny = ny
     patch_data%nz = nz
     patch_data%px = px
     patch_data%py = py
     patch_data%pz = pz
     patch_data%x1 = x1
     patch_data%xn = xn
     patch_data%y1 = y1
     patch_data%yn = yn
     patch_data%z1 = z1
     patch_data%zn = zn
     patch_data%bx1 = TRIM(bx1)
     patch_data%bxn = TRIM(bxn)
     patch_data%by1 = TRIM(by1)
     patch_data%byn = TRIM(byn)
     patch_data%bz1 = TRIM(bz1)
     patch_data%bzn = TRIM(bzn)
     patch_data%simtime = simtime
!     patch_data%dt = dt
!     patch_data%stability = '  startup '
!     patch_data%ijkstab = [0,0,0]  
     patch_data%ax = patch_data%nx/patch_data%px
     patch_data%ay = patch_data%ny/patch_data%py
     patch_data%az = patch_data%nz/patch_data%pz   
     patch_data%dx = (patch_data%xn-patch_data%x1)/REAL(patch_data%nx,KIND=c_double)
     patch_data%dy = (patch_data%yn-patch_data%y1)/REAL(patch_data%ny,KIND=c_double)    
     patch_data%dz = (patch_data%zn-patch_data%z1)/REAL(patch_data%nz,KIND=c_double) 
!     patch_data%isymX = 1
!     patch_data%isymY = 1
!     patch_data%isymZ = 1
     IF ( (patch_data%bx1 == 'SYMM') .OR. (patch_data%bxn == 'SYMM') ) patch_data%isymX = -1
     IF ( (patch_data%by1 == 'SYMM') .OR. (patch_data%byn == 'SYMM') ) patch_data%isymY = -1
     IF ( (patch_data%bz1 == 'SYMM') .OR. (patch_data%bzn == 'SYMM') ) patch_data%isymZ = -1
     SELECT CASE(patch_data%bx1)
     CASE('PERI')
      patch_data%periodicx = .TRUE.
     CASE DEFAULT
      patch_data%periodicx = .FALSE.
     END SELECT
     SELECT CASE(patch_data%by1)
     CASE('PERI')
      patch_data%periodicy = .TRUE.
     CASE DEFAULT
      patch_data%periodicy = .FALSE.
     END SELECT
     SELECT CASE(patch_data%bz1)
     CASE('PERI')
      patch_data%periodicz = .TRUE.
     CASE DEFAULT
      patch_data%periodicz = .FALSE.
   END SELECT

   
   END SUBROUTINE setup_patch

   SUBROUTINE remove_patch(patch_data)
    CLASS(patch_type), INTENT(OUT) :: patch_data  ! auto deallocation of all components on entry
    IF (.FALSE.) patch_data%coordsys = 0          ! to keep the compiler from warning about intent(out)
   END SUBROUTINE remove_patch
      
!===================================================================================================
 END MODULE LES_patch
!===================================================================================================
