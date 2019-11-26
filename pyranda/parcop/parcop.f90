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
MODULE parcop

  USE LES_objects
  USE LES_comm, ONLY : LES_comm_world
  USE LES_compact_operators, ONLY : d1x,d1y,d1z,d4x,d4y,d4z,d8x,d8y,d8z
  USE LES_operators, ONLY : div,grad,Laplacian
  USE LES_operators, ONLY : curl,cross,filter,ring,ringV,filterGdir
  USE LES_operators, ONLY : get_rands_normal, filtRands
  USE LES_explicit
  USE LES_ghost
  USE LES_timers
  
  LOGICAL :: use_explicit = .True.    ! Use new exp. der routines?
  LOGICAL :: ghost_explicit = .False.  ! Do ghost comm for exp. routines?
  
  CONTAINS


    SUBROUTINE setup(patch,level,COMM, &
         & nx,ny,nz,px,py,pz,coordsys, &
         & x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: patch,level
      INTEGER,               INTENT(IN) :: COMM
      INTEGER,               INTENT(IN) :: nx,ny,nz,px,py,pz
      INTEGER,               INTENT(IN) :: coordsys
      REAL(kind=8),   INTENT(IN)   :: x1,xn,y1,yn,z1,zn
      CHARACTER(LEN=*), INTENT(IN) :: bx1,bxn,by1,byn,bz1,bzn
      
      REAL(kind=8)                 :: x1f,xnf,y1f,ynf,z1f,znf
      INTEGER                      :: color,key
      REAL(kind=8)                 :: simtime
      REAL(kind=8)                 :: dx,dy,dz
      

      color = 0
      key = 0
      simtime = 0.0D0

      LES_comm_world = COMM

      ! Convert to faces
      dx = (xn-x1) / DBLE( MAX(nx-1,1) )
      dy = (yn-y1) / DBLE( MAX(ny-1,1) )
      dz = (zn-z1) / DBLE( MAX(nz-1,1) )

      x1f = x1 - dx/2.0D0
      xnf = xn + dx/2.0D0
      y1f = y1 - dy/2.0D0
      ynf = yn + dy/2.0D0
      z1f = z1 - dz/2.0D0
      znf = zn + dz/2.0D0
      
      CALL setup_objects(patch,level,color,key,coordsys,nx,ny,nz,px,py,pz,x1f,xnf,y1f,ynf,z1f,znf,bx1,bxn,by1,byn,bz1,bzn,simtime)
      
      
    END SUBROUTINE setup


    SUBROUTINE setup_mesh(patch,level)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: patch,level

      CALL setup_mesh_data(patch,level)
      
    END SUBROUTINE setup_mesh


    SUBROUTINE setup_mesh_x3(patch,level,x1,x2,x3,meshPer) !mesh_perX,mesh_perY,mesh_perZ)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: patch,level
      REAL(kind=8), DIMENSION(:,:,:), INTENT(IN) :: x1,x2,x3
      LOGICAL, INTENT(IN) :: meshPer !mesh_perX,mesh_perY,mesh_perZ
      
      CALL setup_mesh_data_x3(patch,level,x1,x2,x3,meshPer) !mesh_perX,mesh_perY,mesh_perZ)
      
    END SUBROUTINE setup_mesh_x3


    SUBROUTINE getVar(dxg,vname,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      CHARACTER(LEN=3), INTENT(IN) :: vname
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg

      SELECT CASE(TRIM(vname))
      CASE('x')
         dxg = mesh_ptr%xgrid
      CASE('y')
         dxg = mesh_ptr%ygrid
      CASE('z')
         dxg = mesh_ptr%zgrid
      CASE('d1')
         dxg = mesh_ptr%d1
      CASE('d2')
         dxg = mesh_ptr%d2
      CASE('d3')
         dxg = mesh_ptr%d3
      CASE('dAx')
         dxg = mesh_ptr%dAdx
      CASE('dAy')
         dxg = mesh_ptr%dAdy
      CASE('dAz')
         dxg = mesh_ptr%dAdz
      CASE('dBx')
         dxg = mesh_ptr%dBdx
      CASE('dBy')
         dxg = mesh_ptr%dBdy
      CASE('dBz')
         dxg = mesh_ptr%dBdz
      CASE('dCx')
         dxg = mesh_ptr%dCdx
      CASE('dCy')
         dxg = mesh_ptr%dCdy
      CASE('dCz')
         dxg = mesh_ptr%dCdz
      CASE('dtJ')
         dxg = mesh_ptr%detxyz
      CASE DEFAULT
         print*,"You requested variable: ", TRIM(vname)
         print*,".... I cant find that"
      END SELECT
                          
    END SUBROUTINE getVar
    
    SUBROUTINE xGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg
      
      dxg = mesh_ptr%xgrid      
    END SUBROUTINE xGrid

    SUBROUTINE yGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg      

      dxg = mesh_ptr%ygrid      
    END SUBROUTINE yGrid

    SUBROUTINE zGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg      

      dxg = mesh_ptr%zgrid      
    END SUBROUTINE zGrid

        SUBROUTINE dxGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg
      
      dxg = mesh_ptr%d1     
    END SUBROUTINE dxGrid

    SUBROUTINE dyGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg      

      dxg = mesh_ptr%d2     
    END SUBROUTINE dyGrid

    SUBROUTINE dzGrid(dxg,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dxg      

      dxg = mesh_ptr%d3      
    END SUBROUTINE dzGrid


    SUBROUTINE mesh_getCellVol(CellVol,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      REAL(KIND=8), DIMENSION(nx,ny,nz),INTENT(out) :: CellVol
      CellVol = mesh_ptr%CellVol            
    END SUBROUTINE mesh_getCellVol

    SUBROUTINE mesh_getGridLen(GridLen,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      REAL(KIND=8), DIMENSION(nx,ny,nz),INTENT(out) :: GridLen
      GridLen = mesh_ptr%GridLen
    END SUBROUTINE mesh_getGridLen


    

    SUBROUTINE set_patch(patch,level)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: patch,level
      CALL point_to_objects(patch,level)
    END SUBROUTINE set_patch

    SUBROUTINE divergence(fx,fy,fz,val,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,    INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(in)  :: fx,fy,fz 
      real(kind=8), dimension(nx,ny,nz),intent(out) :: val

      CALL div(fx,fy,fz,val)
      
      
    END SUBROUTINE divergence
    
    SUBROUTINE ddx(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      real(kind=8), dimension(-2:nx+3,ny,nz) :: gval,gdval
      
      ! Stanfard 10th order derivative
      if ( .not. use_explicit ) then
         CALL d1x(val,dval)
      else

         if (ghost_explicit) then
            ! Bin for exp. data      
            gval(1:nx,:,:) = val
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghostx(1,gval)            
            CALL der1e6(gval,gdval,compact_ptr%dx,1,0,0)        
            dval = gdval(1:nx,:,:)
         else         
            CALL der1e6(val,dval,compact_ptr%dx,1,0,0)         
         end if                     
      end if
      
    END SUBROUTINE ddx


    SUBROUTINE ddx2nd(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      

      ! Stanfard 10th order derivative
      !CALL d1x(val,dval)

      ! Bin for exp. data
      real(kind=8), dimension(0:nx+1,ny,nz) :: gval,gdval
      gval(1:nx,:,:) = val

      
      ! Explicit derivative (6th order)
      ! option to handle ghost data
      if (ghost_explicit) then
         CALL ghostx(1,gval)
      end if
      CALL der1e2(gval,gdval,compact_ptr%dx,1,0,0)
      
      dval = gdval(1:nx,:,:) 
      
      
    END SUBROUTINE ddx2nd

    
    SUBROUTINE ddy(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      real(kind=8), dimension(nx,-2:ny+3,nz) :: gval,gdval
      
      if ( .not. use_explicit ) then
         CALL d1y(val,dval)
      else

         if (ghost_explicit) then
            ! Bin for exp. data
            gval(:,1:ny,:) = val
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data
            
            CALL ghosty(1,gval)
            CALL der1e6(gval,gdval,compact_ptr%dy,2,0,0)     
            dval = gdval(:,1:ny,:)
         else
            CALL der1e6(val,dval,compact_ptr%dy,2,0,0)
         end if
         
      end if
      
    END SUBROUTINE ddy
    
    SUBROUTINE ddz(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      real(kind=8), dimension(nx,ny,-2:nz+3) :: gval,gdval
      
      if ( .not. use_explicit ) then
         CALL d1z(val,dval)
      else
         if (ghost_explicit) then
            ! Bin for exp. data

            gval(:,:,1:nz) = val
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghostz(1,gval)
            CALL der1e6(gval,gdval,compact_ptr%dz,3,0,0)
            dval = gdval(:,:,1:nz)
         else
            CALL der1e6(val,dval,compact_ptr%dz,3,0,0)
         end if
         
      end if
      
    END SUBROUTINE ddz

    SUBROUTINE dd4x(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      real(kind=8), dimension(-2:nx+3,ny,nz) :: gval,gdval
      
      if ( .not. use_explicit ) then
         ! This uses compact.f
         !CALL d4x(val,dval)
         CALL d8x(val,dval)
      else

         if (ghost_explicit) then
            ! New explicit.f method
            ! Bin for exp. data
            gval(1:nx,:,:) = val
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghostx(1,gval)
            
            CALL der4e4(gval,gdval,compact_ptr%dx,1,0,0)
            dval = gdval(1:nx,:,:)
         else
            CALL der4e4(val,dval,compact_ptr%dx,1,0,0)
         end if
         
      end if
         
    END SUBROUTINE dd4x
      
            
    SUBROUTINE dd4y(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      real(kind=8), dimension(nx,-2:ny+3,nz) :: gval,gdval
      
      if ( .not. use_explicit ) then
         !CALL d4y(val,dval)
         CALL d8y(val,dval)
      else
         if (ghost_explicit) then
            ! New explicit.f method
            ! Bin for exp. data
            gval(:,1:ny,:) = val
         
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghosty(1,gval)
            CALL der4e4(gval,gdval,compact_ptr%dy,2,0,0)        
            dval = gdval(:,1:ny,:)
         else
            CALL der4e4(val,dval,compact_ptr%dy,2,0,0)
         end if
      end if
      
      
    END SUBROUTINE dd4y

    SUBROUTINE dd4z(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      real(kind=8), dimension(nx,ny,-2:nz+3) :: gval,gdval
      
      if ( .not. use_explicit ) then
         !CALL d4z(val,dval)
         CALL d8z(val,dval)
      else

         if (ghost_explicit) then
            ! New explicit.f method
            ! Bin for exp. data
            gval(:,:,1:nz) = val         
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghostz(1,gval)
            CALL der4e4(gval,gdval,compact_ptr%dz,3,0,0)      
            dval = gdval(:,:,1:nz)
         else
            CALL der4e4(val,dval,compact_ptr%dz,3,0,0)
         end if
      end if
      
    END SUBROUTINE dd4z

    

    SUBROUTINE plaplacian(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      
      CALL Laplacian(val,dval)

    END SUBROUTINE plaplacian
   
    SUBROUTINE pRing(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      
      dval = ring( val )

    END SUBROUTINE pRing

    SUBROUTINE pRingV(vx,vy,vz,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: vx,vy,vz
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      
      dval = ringV( vx, vy, vz )

    END SUBROUTINE pRingV
 
    

    SUBROUTINE sFilter(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      real(kind=8), dimension(nx,ny,nz) :: tmp
      CHARACTER(LEN=8), PARAMETER :: filtype='spectral'
      real(kind=8), dimension(-2:nx+3,-2:ny+3,-2:nz+3) :: gval,gdval

      !$omp target data map(alloc:tmp) 
      
      ! This uses compact.f90
      if ( .not. use_explicit ) then
         CALL filter(filtype,val,dval)
      else

         if (ghost_explicit) then
            ! Bin for exp. data
            gval(1:nx,1:ny,1:nz) = val
            
            ! Explicit derivative (6th order)
            ! option to handle ghost data
            
            CALL ghostx(1,gval)
            CALL ghosty(1,gval)
            CALL ghostz(1,gval)
            
            
            CALL filte6(val,dval,1,0,0,nx,ny,nz)
            CALL filte6(dval,tmp,2,0,0,nx,ny,nz)
            CALL filte6(tmp,dval,3,0,0,nx,ny,nz)
         
            dval = gdval(1:nx,1:ny,1:nz)
         else
            
            CALL filte6(val,dval,1,0,0,nx,ny,nz)
            CALL filte6(dval,tmp,2,0,0,nx,ny,nz)
            CALL filte6(tmp,dval,3,0,0,nx,ny,nz)
         endif
      end if

   !$omp end target data
      
    END SUBROUTINE sFilter

    SUBROUTINE gFilter(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      real(kind=8), dimension(nx,ny,nz) :: tmp
      CHARACTER(LEN=6), PARAMETER :: filtype='smooth'
      real(kind=8), dimension(-2:nx+3,-2:ny+3,-2:nz+3) :: gval,gdval

      !$omp target data map(alloc:tmp)      
      
      if ( .not. use_explicit ) then
         ! Compact version
         CALL filter(filtype,val,dval)
      else

         if (ghost_explicit) then
            ! Bin for exp. data
            gval(1:nx,1:ny,1:nz) = val
      
            ! Explicit derivative (6th order)
            ! option to handle ghost data

            CALL ghostx(1,gval)
            CALL ghosty(1,gval)
            CALL ghostz(1,gval)
            
            CALL gfilt3(val,dval,1,0,0,nx,ny,nz)
            CALL gfilt3(dval,tmp,2,0,0,nx,ny,nz)
            CALL gfilt3(tmp,dval,3,0,0,nx,ny,nz)         
            
            dval = gdval(1:nx,:,:)
         else
            
            CALL gfilt3(val,dval,1,0,0,nx,ny,nz)
            CALL gfilt3(dval,tmp,2,0,0,nx,ny,nz)
            CALL gfilt3(tmp,dval,3,0,0,nx,ny,nz)         
         end if
         
      end if

      !$omp end target data
      
    END SUBROUTINE gFilter
    
    SUBROUTINE gFilterDir(val,dval,nx,ny,nz,dir)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz,dir
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      CHARACTER(LEN=6), PARAMETER :: filtype='smooth'

      CALL filterGdir(filtype,val,dval,dir)

    END SUBROUTINE gFilterDir


    SUBROUTINE gradS(val,val1,val2,val3,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz),intent(in) :: val
      real(kind=8), dimension(nx,ny,nz), intent(out) :: val1,val2,val3
    
      CALL grad(val,val1,val2,val3)

    END SUBROUTINE gradS


   

 !! Communication objects from Python
    SUBROUTINE commFromPy( patch, level, comm_list )
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: patch,level
      INTEGER, DIMENSION(7), INTENT(IN) :: comm_list

      CALL point_to_objects(patch,level)
      
      comm_ptr%xyzcom = comm_list(1)
      comm_ptr%xycom  = comm_list(2)
      comm_ptr%xzcom  = comm_list(3)
      comm_ptr%yzcom  = comm_list(4)
      comm_ptr%xcom   = comm_list(5)
      comm_ptr%ycom   = comm_list(6)
      comm_ptr%zcom   = comm_list(7)
      
      ! Need to add _hi _lo comms for ghost data

      
    END SUBROUTINE commFromPy


    SUBROUTINE commx(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%xcom
    END SUBROUTINE commx

    SUBROUTINE commy(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%ycom
    END SUBROUTINE commy

    SUBROUTINE commz(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%zcom
    END SUBROUTINE commz

    SUBROUTINE commxy(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%xycom
    END SUBROUTINE commxy

    SUBROUTINE commxz(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%xzcom
    END SUBROUTINE commxz

    SUBROUTINE commyz(COMM)
      IMPLICIT NONE
      INTEGER,               INTENT(OUT) :: COMM      
      COMM = comm_ptr%yzcom
    END SUBROUTINE commyz



  SUBROUTINE TBL_get_rands(rands,ny,nz,Nbuff,time_seed)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(4,ny+Nbuff,nz+Nbuff),INTENT(OUT) :: rands
    INTEGER, INTENT(IN) :: ny,nz,Nbuff
    INTEGER, INTENT(IN) :: time_seed 

    CALL get_rands_normal( rands, ny, nz, Nbuff, time_seed)

  END SUBROUTINE TBL_get_rands


  SUBROUTINE TBL_filter(Nspan,Ni,No,Nbuff, &
       bmnI,bmnO,rands,vfilt, &
       ny, nz, ay, az, iy1, iz1,y_r )        
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nspan,Ni,No,Nbuff
    INTEGER, INTENT(IN) :: ny,nz,iy1,iz1,ay,az
    DOUBLE PRECISION, INTENT(IN) :: y_r(ay)
    DOUBLE PRECISION, INTENT(IN) ::  rands(ny+Nbuff,nz+Nbuff)
    DOUBLE PRECISION, INTENT(IN) ::  bmnI(2*Ni+1,2*Nspan+1), bmnO(2*No+1,2*Nspan+1)
    DOUBLE PRECISION, DIMENSION(ay,az), INTENT(OUT) :: vfilt(ay,az)
    INTEGER :: N1,N2
    INTEGER :: j,k,m,n,mm,nn
    INTEGER :: mF,mG,nF,nG

    CALL filtRands( Nspan,Ni,No,Nbuff, &
         bmnI,bmnO,rands,vfilt, &
         ny, nz, ay, az, iy1, iz1,y_r )
    
  END SUBROUTINE TBL_filter


END MODULE parcop

