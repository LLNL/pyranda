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
 MODULE LES_mesh ! MESH VARIABLES
!===================================================================================================  
  USE iso_c_binding
  USE LES_patch, ONLY : patch_type
  USE LES_comm, ONLY : comm_type
  USE LES_compact, ONLY : compact_type,compact_op1,mesh1_type,comm1_type,control_data,compact_weight_d1,MPI_PROC_NULL,compact_op1_d1
  !USE LES_cb_dd, ONLY : lua_cb_dd
  !#ifdef LUA
  !USE wkt_mesh
  ! #endif
  IMPLICIT NONE

  logical, parameter :: custom=.false.

  TYPE mesh_type ! fully allocatable; scalars are allocated on assignment
    LOGICAL(c_bool) :: x1proc,y1proc,z1proc,xnproc,ynproc,znproc      ! Process contains boundary?
    INTEGER(c_int) :: num_zones                                       ! Number of local zones
    INTEGER(c_int),  DIMENSION(:)    , ALLOCATABLE :: zone_index      ! 1D zonal index
    INTEGER(c_int),  DIMENSION(:,:,:), ALLOCATABLE :: refine          ! Refinement flag for SAMRAI
    INTEGER(c_int),  DIMENSION(:)    , ALLOCATABLE :: ix              ! Global x index of grid 
    INTEGER(c_int),  DIMENSION(:)    , ALLOCATABLE :: iy              ! Global y index of grid
    INTEGER(c_int),  DIMENSION(:)    , ALLOCATABLE :: iz              ! Global z index of grid
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: xgrid           ! x or radial location (radius => xgrid)
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: ygrid           ! y or theta location  (theta  => ygrid)
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: zgrid           ! z or phi location    (phi    => zgrid)
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: d1              ! Dimensional grid spacing in x or radius
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: d2              ! Dimensional grid spacing in y or theta 
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: d3              ! Dimensional grid spacing in z or phi 
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: GridLen         ! Minimum grid spacing
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: CellVol         ! Cell volume
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: CellVolS        ! Sharp-filtered cell volume
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: CellVolG        ! Gaussian-filtered cell volume
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: detxyz          ! Determinant of Jacobian
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: dAdx,dAdy,dAdz  ! A-row of inverse Jacobian tensor
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: dBdx,dBdy,dBdz  ! B-row of inverse Jacobian tensor
    REAL(c_double),  DIMENSION(:,:,:), ALLOCATABLE :: dCdx,dCdy,dCdz  ! C-row of inverse Jacobian tensor
    CONTAINS
     PROCEDURE :: setup => setup_mesh
     PROCEDURE :: remove => remove_mesh
!     FINAL :: finalize_mesh  ! not needed here
  END TYPE mesh_type
  
!===================================================================================================
  CONTAINS
!===================================================================================================

   SUBROUTINE setup_mesh(mesh_data,patch_data,comm_data,compact_data,xpy,ypy,zpy,custom_per) !custom_periodicX,custom_periodicY,custom_periodicZ)
    IMPLICIT NONE
    CLASS(mesh_type),  INTENT(OUT) :: mesh_data ! Auto deallocation of all components on entry
    CLASS(patch_type), INTENT(IN)  :: patch_data
    CLASS(comm_type),  INTENT(IN)  :: comm_data
    CLASS(compact_type),  INTENT(IN)  :: compact_data
    !REAL(c_double), DIMENSION(patch_data%nx,patch_data%ny,patch_data%nz), INTENT(IN), OPTIONAL :: xpy,ypy,zpy
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: xpy
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: ypy
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: zpy
    logical,INTENT(IN),OPTIONAL  :: custom_per
    !logical,INTENT(IN),OPTIONAL  :: custom_periodicX
    !logical,INTENT(IN),OPTIONAL  :: custom_periodicY
    !logical,INTENT(IN),OPTIONAL  :: custom_periodicZ
    !TYPE(compact_op1)    :: custom_op
    TYPE(compact_op1_d1) :: custom_d1
    TYPE(mesh1_type) :: custom_mesh
    TYPE(comm1_type) :: custom_comm
    integer :: custom_lo,custom_hi
    logical :: custom_periodic
    logical(c_bool) :: custom_nullop
    REAL(c_double), PARAMETER :: half=0.5_c_double, pi=acos(-1.0_c_double)
    REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: dxdA,dxdB,dxdC	!  X-row of Jacobian tensor
    REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: dydA,dydB,dydC	!  Y-row of Jacobian tensor
    REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: dzdA,dzdB,dzdC	!  Z-row of Jacobian tensor
    REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: tmp                !  work array
    REAL(c_double) :: dA,dB,dC
    LOGICAL :: tmpPx,tmpPy,tmpPz
    INTEGER(c_int) :: i,j,k,flag,filnum,xasym,yasym,zasym
        
     IF ( .NOT. ALLOCATED(mesh_data%zone_index) ) THEN 
     mesh_data%num_zones = patch_data%ax*patch_data%ay*patch_data%az 
     ALLOCATE(mesh_data%zone_index(mesh_data%num_zones)); FORALL(i=1:mesh_data%num_zones) mesh_data%zone_index(i) = i-1
     ALLOCATE(mesh_data%refine(patch_data%ax,patch_data%ay,patch_data%az))  
     ALLOCATE(mesh_data%ix(patch_data%ax))
     ALLOCATE(mesh_data%iy(patch_data%ay))
     ALLOCATE(mesh_data%iz(patch_data%az))
     ALLOCATE(mesh_data%xgrid(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%ygrid(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%zgrid(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%d1(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%d2(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%d3(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%GridLen(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%CellVol(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%CellVolS(patch_data%ax,patch_data%ay,patch_data%az))
     ALLOCATE(mesh_data%CellVolG(patch_data%ax,patch_data%ay,patch_data%az))
     SELECT CASE(patch_data%coordsys)
     CASE(3)
       ALLOCATE(mesh_data%detxyz(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dAdx(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dAdy(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dAdz(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dBdx(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dBdy(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dBdz(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dCdx(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dCdy(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(mesh_data%dCdz(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dxdA(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dxdB(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dxdC(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dydA(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dydB(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dydC(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dzdA(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dzdB(patch_data%ax,patch_data%ay,patch_data%az))
       ALLOCATE(dzdC(patch_data%ax,patch_data%ay,patch_data%az))
     END SELECT
     ENDIF
     ASSOCIATE(coordsys  => patch_data%coordsys,  simtime   => patch_data%simtime,                                      &
               nx        => patch_data%nx,        ny        => patch_data%ny,        nz        => patch_data%nz,        &
               px        => patch_data%px,        py        => patch_data%py,        pz        => patch_data%pz,        &
               x1        => patch_data%x1,        y1        => patch_data%y1,        z1        => patch_data%z1,        &
               xn        => patch_data%xn,        yn        => patch_data%yn,        zn        => patch_data%zn,        &
               bx1       => patch_data%bx1,       by1       => patch_data%by1,       bz1       => patch_data%bz1,       &
               bxn       => patch_data%bxn,       byn       => patch_data%byn,       bzn       => patch_data%bzn,       &
               ax        => patch_data%ax,        ay        => patch_data%ay,        az        => patch_data%az,        &
               dx        => patch_data%dx,        dy        => patch_data%dy,        dz        => patch_data%dz,        &
               periodicx => patch_data%periodicx, periodicy => patch_data%periodicy, periodicz => patch_data%periodicz, &
               xcom_id   => comm_data%xcom_id,    ycom_id   => comm_data%ycom_id,    zcom_id   => comm_data%zcom_id,    &
               x1proc    => mesh_data%x1proc,     y1proc    => mesh_data%y1proc,     z1proc    => mesh_data%z1proc,     &
               xnproc    => mesh_data%xnproc,     ynproc    => mesh_data%ynproc,     znproc    => mesh_data%znproc,     &
               ix        => mesh_data%ix,         iy        => mesh_data%iy,         iz        => mesh_data%iz,         &
               xgrid     => mesh_data%xgrid,      ygrid     => mesh_data%ygrid,      zgrid     => mesh_data%zgrid,      &
               d1        => mesh_data%d1,         d2        => mesh_data%d2,         d3        => mesh_data%d3,         &
               GridLen   => mesh_data%GridLen,    CellVol   => mesh_data%CellVol,    CellVolS  => mesh_data%CellVolS,   &
	       CellVolG  => mesh_data%CellVolG,   detxyz    => mesh_data%detxyz,                                        &
               dAdx      => mesh_data%dAdx,       dAdy      => mesh_data%dAdy,       dAdz      => mesh_data%dAdz,       &
               dBdx      => mesh_data%dBdx,       dBdy      => mesh_data%dBdy,       dBdz      => mesh_data%dBdz,       &
               dCdx      => mesh_data%dCdx,       dCdy      => mesh_data%dCdy,       dCdz      => mesh_data%dCdz,       &
	       radius    => mesh_data%xgrid,      theta     => mesh_data%ygrid)
     ix = xcom_id*ax + [ (i,i=1,ax) ]
     iy = ycom_id*ay + [ (j,j=1,ay) ]
     iz = zcom_id*az + [ (k,k=1,az) ]
     x1proc = ( ix(1)  == 1  .and. (.not. periodicx) .and. nx > 1 )
     y1proc = ( iy(1)  == 1  .and. (.not. periodicy) .and. ny > 1 )
     z1proc = ( iz(1)  == 1  .and. (.not. periodicz) .and. nz > 1 )
     xnproc = ( ix(ax) == nx .and. (.not. periodicx) .and. nx > 1 )
     ynproc = ( iy(ay) == ny .and. (.not. periodicy) .and. ny > 1 )
     znproc = ( iz(az) == nz .and. (.not. periodicz) .and. nz > 1 )

     IF ( PRESENT(xpy) ) THEN  ! All or nothing on passing grid explicitly
        xgrid = xpy
        ygrid = ypy
        zgrid = zpy
     ELSE
        FORALL (i=1:ax,j=1:ay,k=1:az)
           xgrid(i,j,k) = x1+REAL(2*ix(i)-1,KIND=c_double)*half*dx ! AKA radius
           ygrid(i,j,k) = y1+REAL(2*iy(j)-1,KIND=c_double)*half*dy ! AKA theta
           zgrid(i,j,k) = z1+REAL(2*iz(k)-1,KIND=c_double)*half*dz ! AKA phi
        END FORALL
     END IF

     
     
! CELL VOLUMES AND MINIMUM GRID SPACINGS -----------------------------------------------------------     
     SELECT CASE(coordsys)
     CASE(0) ! Cartesian
       CellVol = dx*dy*dz 
       SELECT CASE(nx)
       CASE(1)
         d1 = MAX(dy,dz)
       CASE DEFAULT	 
         d1 = dx
       END SELECT
       SELECT CASE(ny)
       CASE(1)
         d2 = MAX(dx,dz)
       CASE DEFAULT	 
         d2 = dy
       END SELECT
       SELECT CASE(nz)
       CASE(1)
         d3 = MAX(dx,dy)
       CASE DEFAULT	 
         d3 = dz
       END SELECT  	 
       GridLen = MIN(d1,d2,d3)
       RETURN      
     CASE(1) ! Cylindrical
       d1 =        dx
       d2 = radius*dy
       d3 =        dz 
       SELECT CASE(ny)
       CASE(1)
         CellVol = d1*2.0D0*pi*radius*d3
	 GridLen = MIN(d1,d3)
       CASE DEFAULT
         CellVol = d1*d2*d3
         GridLen = MIN(d1,d2,d3)  
       END SELECT
     CASE(2) ! Spherical
       d1 =                   dx
       d2 =            radius*dy
       d3 = radius*SIN(theta)*dz    
       CellVol = d1*d2*d3
       GridLen = MIN(d1,d2,d3)  
       SELECT CASE(ny)
       CASE(1)
         SELECT CASE(nz)
         CASE(1)
           CellVol = d1*4.0D0*pi*radius**2
           GridLen = d1
	 END SELECT ! nz
       END SELECT ! ny
     CASE(3) ! Curvilinear
     
       dA = dx
       dB = dy
       dC = dz
       
       
     if( custom_per ) then ! custom 1st derivative with unknown bcs

       custom_periodic = .false. ; custom_nullop = .false.
       custom_lo = comm_data%xcom_lo ; custom_hi = comm_data%xcom_hi
       if( comm_data%xcom_id == 0 ) custom_lo = MPI_PROC_NULL
       if( comm_data%xcom_id == comm_data%xcom_np-1 ) custom_hi = MPI_PROC_NULL
       custom_comm = comm1_type(custom_periodic,comm_data%xcom,comm_data%xcom_np,comm_data%xcom_id, &
         custom_lo,custom_hi,comm_data%xrange)
!       custom_mesh = mesh1_type(coordsys,ax,nx,dx,x1,xn,'NONE','NONE')
       custom_mesh = mesh1_type(ax,nx,dx,x1,xn,'NONE','NONE')
       CALL custom_d1%setup(compact_weight_d1(control_data%d1spec),custom_comm,custom_mesh,[0,0],custom_nullop)


       CALL custom_d1%evalx(ax,ay,az,xgrid,dxdA) 
       CALL custom_d1%evalx(ax,ay,az,ygrid,dydA) 
       CALL custom_d1%evalx(ax,ay,az,zgrid,dzdA) 
       dxdA = dxdA / dA
       dydA = dydA / dA
       dzdA = dzdA / dA

       custom_lo = comm_data%ycom_lo ; custom_hi = comm_data%ycom_hi
       if( comm_data%ycom_id == 0 ) custom_lo = MPI_PROC_NULL
       if( comm_data%ycom_id == comm_data%ycom_np-1 ) custom_hi = MPI_PROC_NULL
       custom_comm = comm1_type(custom_periodic,comm_data%ycom,comm_data%ycom_np,comm_data%ycom_id, &
         custom_lo,custom_hi,comm_data%yrange)
!       custom_mesh = mesh1_type(coordsys,ay,ny,dy,y1,yn,'NONE','NONE')
       custom_mesh = mesh1_type(ay,ny,dy,y1,yn,'NONE','NONE')
       CALL custom_d1%setup(compact_weight_d1(control_data%d1spec),custom_comm,custom_mesh,[0,0],custom_nullop)
       CALL custom_d1%evaly(ax,ay,az,xgrid,dxdB)
       CALL custom_d1%evaly(ax,ay,az,ygrid,dydB)
       CALL custom_d1%evaly(ax,ay,az,zgrid,dzdB)
       
       dxdB = dxdB / dB
       dydB = dydB / dB
       dzdB = dzdB / dB


       custom_lo = comm_data%zcom_lo ; custom_hi = comm_data%zcom_hi
       if( comm_data%zcom_id == 0 ) custom_lo = MPI_PROC_NULL
       if( comm_data%zcom_id == comm_data%zcom_np-1 ) custom_hi = MPI_PROC_NULL
       custom_comm = comm1_type(custom_periodic,comm_data%zcom,comm_data%zcom_np,comm_data%zcom_id, &
         custom_lo,custom_hi,comm_data%zrange)
!       custom_mesh = mesh1_type(coordsys,az,nz,dz,z1,zn,'NONE','NONE')
       custom_mesh = mesh1_type(az,nz,dz,z1,zn,'NONE','NONE')
       CALL custom_d1%setup(compact_weight_d1(control_data%d1spec),custom_comm,custom_mesh,[0,0],custom_nullop)
       CALL custom_d1%evalz(ax,ay,az,xgrid,dxdC)
       CALL custom_d1%evalz(ax,ay,az,ygrid,dydC)
       CALL custom_d1%evalz(ax,ay,az,zgrid,dzdC)
       
       dxdC = dxdC / dC
       dydC = dydC / dC
       dzdC = dzdC / dC

       call custom_d1%remove()

     else ! use stock routines for scalar

       CALL compact_data%d1x(1)%evalx(ax,ay,az,xgrid,dxdA) 
       CALL compact_data%d1x(1)%evalx(ax,ay,az,ygrid,dydA) 
       CALL compact_data%d1x(1)%evalx(ax,ay,az,zgrid,dzdA) 
       dxdA = dxdA / dA
       dydA = dydA / dA
       dzdA = dzdA / dA

       CALL compact_data%d1y(1)%evaly(ax,ay,az,xgrid,dxdB)
       CALL compact_data%d1y(1)%evaly(ax,ay,az,ygrid,dydB)
       CALL compact_data%d1y(1)%evaly(ax,ay,az,zgrid,dzdB)       
       dxdB = dxdB / dB
       dydB = dydB / dB
       dzdB = dzdB / dB

       CALL compact_data%d1z(1)%evalz(ax,ay,az,xgrid,dxdC)
       CALL compact_data%d1z(1)%evalz(ax,ay,az,ygrid,dydC)
       CALL compact_data%d1z(1)%evalz(ax,ay,az,zgrid,dzdC)       
       dxdC = dxdC / dC
       dydC = dydC / dC
       dzdC = dzdC / dC

     endif ! custom

       !  Special cases for 2d grids
       IF (nx .EQ. 1) dxdA = 1.0d0
       IF (ny .EQ. 1) dydB = 1.0d0
       IF (nz .EQ. 1) dzdC = 1.0d0
     
       !  Get the determinant of J
       detxyz = -dxdC*dydB*dzdA+dxdB*dydC*dzdA+dxdC*dydA*dzdB &
                -dxdA*dydC*dzdB-dxdB*dydA*dzdC+dxdA*dydB*dzdC
     
       !  Get inv(J)
       dAdx = (-dydC*dzdB + dydB*dzdC)/detxyz
       dAdy = ( dxdC*dzdB - dxdB*dzdC)/detxyz
       dAdz = (-dxdC*dydB + dxdB*dydC)/detxyz
     
       dBdx = ( dydC*dzdA - dydA*dzdC)/detxyz
       dBdy = (-dxdC*dzdA + dxdA*dzdC)/detxyz
       dBdz = ( dxdC*dydA - dxdA*dydC)/detxyz
     
       dCdx = (-dydB*dzdA + dydA*dzdB)/detxyz
       dCdy = ( dxdB*dzdA - dxdA*dzdB)/detxyz
       dCdz = (-dxdB*dydA + dxdA*dydB)/detxyz

       !  Physical local grid length scale in A,B,and C grid directions
       d1 = sqrt( (dxdA*dA)**2 + (dydA*dA)**2 + (dzdA*dA)**2 ) 
       d2 = sqrt( (dxdB*dB)**2 + (dydB*dB)**2 + (dzdB*dB)**2 )
       d3 = sqrt( (dxdC*dC)**2 + (dydC*dC)**2 + (dzdC*dC)**2 )
       CellVol = d1*d2*d3
       GridLen = MIN(d1,d2,d3)

       DEALLOCATE(dzdC)
       DEALLOCATE(dzdB)
       DEALLOCATE(dzdA)
       DEALLOCATE(dydC)
       DEALLOCATE(dydB)
       DEALLOCATE(dydA)
       DEALLOCATE(dxdC)
       DEALLOCATE(dxdB)
       DEALLOCATE(dxdA)
       
     END SELECT ! coordsys
     
! In order to filter on non-cartesian grids, we need to divide by the filtered cell volume.---------     
     
     SELECT CASE(coordsys)
     CASE(2,3) ! Spherical, Curvilinear for scalars (vectors??)
       ALLOCATE(tmp(patch_data%ax,patch_data%ay,patch_data%az))
        CALL compact_data%sfx(1)%evalx(ax,ay,az,CellVol, CellVolS)
        CALL compact_data%sfy(1)%evaly(ax,ay,az,CellVolS,tmp)
        CALL compact_data%sfz(1)%evalz(ax,ay,az,tmp,     CellVolS)
        CALL compact_data%gfx(1)%evalx(ax,ay,az,CellVol, CellVolG)
        CALL compact_data%gfy(1)%evaly(ax,ay,az,CellVolG,tmp)
        CALL compact_data%gfz(1)%evalz(ax,ay,az,tmp,     CellVolG)
       DEALLOCATE(tmp)
     CASE DEFAULT
       CellVolS = CellVol
       CellVolG = CellVol
     END SELECT
     
     END ASSOCIATE
   END SUBROUTINE setup_mesh

   SUBROUTINE remove_mesh(mesh_data)
    CLASS(mesh_type), INTENT(OUT) :: mesh_data  ! auto deallocation of all components on entry
    IF (.FALSE.) mesh_data%ix(1) = 1            ! to keep the compiler from complaining about intent(out)
   END SUBROUTINE remove_mesh    

!===================================================================================================
 END MODULE LES_mesh
!===================================================================================================
