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
  USE LES_compact_operators, ONLY : d1x,d1y,d1z,d4x,d4y,d4z
  USE LES_operators, ONLY : div,grad,Laplacian
  USE LES_operators, ONLY : curl,cross,filter,ring,ringV,filterGdir
  USE LES_operators, ONLY : get_rands_normal, filtRands
  
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

      CALL d1x(val,dval)

    END SUBROUTINE ddx

    SUBROUTINE ddy(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      
      CALL d1y(val,dval)

    END SUBROUTINE ddy
    
    SUBROUTINE ddz(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      
      CALL d1z(val,dval)

    END SUBROUTINE ddz

    SUBROUTINE dd4x(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      CALL d4x(val,dval)
    END SUBROUTINE dd4x

    SUBROUTINE dd4y(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      CALL d4y(val,dval)
    END SUBROUTINE dd4y

    SUBROUTINE dd4z(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval      
      CALL d4z(val,dval)
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
      CHARACTER(LEN=8), PARAMETER :: filtype='spectral'

      CALL filter(filtype,val,dval)

    END SUBROUTINE sFilter

    SUBROUTINE gFilter(val,dval,nx,ny,nz)
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: nx,ny,nz
      real(kind=8), dimension(nx,ny,nz), intent(in) :: val
      real(kind=8), dimension(nx,ny,nz),intent(out) :: dval
      CHARACTER(LEN=6), PARAMETER :: filtype='smooth'

      CALL filter(filtype,val,dval)

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


  SUBROUTINE flamencoFlux()

    ! Take in pre-ghosted conserved variables and return fluxes on those EOMs
    ! Inputs: iNnx,iNny,iNnz
    !
    !         U (zonal,Conserved Vars) -> rSol3D
    !             size(iNVar,-iNcHalo+1:iNnx+iNcHalo-1,-iNcHalo+1:iNny+iNcHalo-1,-iNcHalo+1:iNnz+iNcHalo-1))
    !         rX,rY,rZ (nodal,grid points)
    !             size(1-iNnHalo:iNnx+iNnHalo,1-iNnHalo:iNny+iNnHalo,1-iNnHalo:iNnz+iNnHalo)
    ! Outputs:
    !         rFlux3D (flux of conserved vars)
    
    implicit none
    integer n
    ! Array sizes, various options
    integer :: iNVar,iNVarComputed,iNcHalo,iNnHalo,iNnx,iNny,iNnz,iLowMachScheme,iNumberOfSpecies,iEquationOfState,iRecons,iInvOrder
    ! Thermodynamic properties
    real(8),allocatable,dimension(:) :: rSpecificGasConstant
    real(8),allocatable,dimension(:,:,:) :: rCv,rCp
    ! X,Y,Z arrays
    real(8),allocatable,dimension(:,:,:) :: rX,rY,rZ
    ! Solution array, flux array
    real(8),allocatable,dimension(:,:,:,:) :: rSol3D,rFlux3D
    
    iNnx = 33 ! Number of nodes in x direction
    iNny = 33 ! Number of nodes in y direction
    iNnz = 33 ! Number of nodes in z direction
    iNumberOfSpecies = 1 ! Number of species
    iNVarComputed = 4+iNumberOfSpecies ! Number of computed variables in solution array
    iNVar = iNVarComputed+4 ! Total number of variables in solution array
    iEquationOfState = 1 ! Equation of state. 1=ideal gas (fixed gamma)
    iRecons = 0 ! Choice of variables to reconstruct at interface. 0=conserved, 1=primitive (preferred for multispecies)
    iInvOrder = 5 ! Reconstruction scheme to use. 1=1st order, 2=minmod, 3=van Leer, 4=superbee, 5=5th order MUSCL
    iLowMachScheme = 1 ! Low Mach correction. 0=off, 1=on
    
    ! Number of halo cells based on reconstruction scheme being used
    select case(iInvOrder)
    case(1)
       iNcHalo=1
    case(2:4)
       iNcHalo=2
    case(5)
       iNcHalo=3
    end select
    ! Number of halo nodes (always 1)
    iNnHalo = 1
    
    ! Array dimensions
    allocate(rX(1-iNnHalo:iNnx+iNnHalo,1-iNnHalo:iNny+iNnHalo,1-iNnHalo:iNnz+iNnHalo))
    allocate(rY(1-iNnHalo:iNnx+iNnHalo,1-iNnHalo:iNny+iNnHalo,1-iNnHalo:iNnz+iNnHalo))
    allocate(rZ(1-iNnHalo:iNnx+iNnHalo,1-iNnHalo:iNny+iNnHalo,1-iNnHalo:iNnz+iNnHalo))
    allocate(rSol3D(iNVar,-iNcHalo+1:iNnx+iNcHalo-1,-iNcHalo+1:iNny+iNcHalo-1,-iNcHalo+1:iNnz+iNcHalo-1))
    allocate(rFlux3D(iNVar,-iNcHalo+1:iNnx+iNcHalo-1,-iNcHalo+1:iNny+iNcHalo-1,-iNcHalo+1:iNnz+iNcHalo-1))
    allocate(rSpecificGasConstant(iNumberOfSpecies))
    allocate(rCv(iNumberOfSpecies,2,5))
    allocate(rCp(iNumberOfSpecies,2,5))
    
    ! (X,Y,Z) coordinates of each node
    rX = 0.d0; rY = 0.d0; rZ = 0.d0
    
    ! Solution array
    rSol3D(1,:,:,:) = 0.d0 ! X momentum (rho*u)
    rSol3D(2,:,:,:) = 0.d0 ! Y momentum (rho*v)
    rSol3D(3,:,:,:) = 0.d0 ! Z momentum (rho*w)
    rSol3D(4,:,:,:) = 0.d0 ! Total energy (rho*E)
    if(iNumberOfSpecies.eq.1) then
       rSol3D(5,:,:,:) = 0.d0 ! Density
    elseif(iNumberOfSpecies.gt.1) then
       do n=1,iNumberOfSpecies
          rSol3D(4+n,:,:,:) = 0.d0 ! Density*Mass Fraction (rho*Y_n)
       end do
       rSol3D(iNVar-3,:,:,:) = 0.d0 ! Density
    end if
    rSol3D(iNVar-2,:,:,:) = 0.d0 ! Pressure
    rSol3D(iNVar-1,:,:,:) = 0.d0 ! Gamma
    rSol3D(iNVar,:,:,:) = 0.d0 ! Temperature
    
    ! Flux array
    rFlux3D = 0.d0 ! Same as rSol3D
    
    ! Thermodynamic properties for fixed-gamma ideal gas
    do n=1,iNumberOfSpecies
       rSpecificGasConstant(n) = 1.d0 ! Specific gas constant for species n
       rCp(n,1,1) = 1.d0 ! Specific heat at constant pressure for species n 
       rCp(n,1,2) = 0.d0; rCp(n,1,3) = 0.d0; rCp(n,1,4) = 0.d0; rCp(n,1,5) = 0.d0 ! For fixed-gamma gas other coefficients are set to zero
       rCp(n,2,1) = 0.d0; rCp(n,2,2) = 0.d0; rCp(n,2,3) = 0.d0; rCp(n,2,4) = 0.d0; rCp(n,2,5) = 0.d0 ! For fixed-gamma gas other coefficients are set to zero
       rCv(n,1,1) = 1.d0 ! Specific heat at constant volume for species n 
       rCv(n,1,2) = 0.d0; rCv(n,1,3) = 0.d0; rCv(n,1,4) = 0.d0; rCv(n,1,5) = 0.d0 ! For fixed-gamma gas other coefficients are set to zero
       rCv(n,2,1) = 0.d0; rCv(n,2,2) = 0.d0; rCv(n,2,3) = 0.d0; rCv(n,2,4) = 0.d0; rCv(n,2,5) = 0.d0 ! For fixed-gamma gas other coefficients are set to zero
       
    end do
    

    call InviscidFlux(rSol3D,rFlux3D,rX,rY,rZ,iNVar,iNVarComputed,iNcHalo,iNnHalo,iNnx,iNny,iNnz,iLowMachScheme,iNumberOfSpecies,iEquationOfState,rSpecificGasConstant,rCv,rCp,iRecons,iInvOrder)
    
    
    
  END SUBROUTINE flamencoFlux

  

END MODULE parcop

