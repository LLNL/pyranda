PROGRAM miniApp
  USE iso_c_binding
  USE LES_comm, ONLY : LES_comm_world
  USE LES_objects
  USE LES_timers
  USE parcop, ONLY : setup,ddx,point_to_objects,setup_mesh,grad,filter,div,ddy,ddz,sfilter,gfilter,dd4x
  USE LES_ghost, ONLY : ghostx,ghosty,ghostz
  IMPLICIT NONE
  INCLUDE "mpif.h"
  
  INTEGER(c_int)                 :: nx,ny,nz,px,py,pz,ax,ay,az,ns
  REAL(c_double)                 :: x1,xn,y1,yn,z1,zn
  CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn
  REAL(c_double)                 :: simtime
  INTEGER(c_int)                 :: world_id,world_np,mpierr
  LOGICAL(c_bool)                :: periodic
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: uxx,uyy,uzz,uxy,uxz,uyz,uyx,uzx,uzy
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: beta,mu,kappa
  REAL(c_double), DIMENSION(:,:,:,:), ALLOCATABLE :: RHS,Y,adiff,RHSp

  REAL(c_double) :: divu,txx,tyy,tzz,txy,txz,tyz
  
  INTEGER :: tt,i,j,k,n,eom
  INTEGER :: t1,t2,clock_rate,clock_max
  DOUBLE PRECISION :: t1c,t2c
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,ii,iterations
  INTEGER :: rank,ierror,procs
  DOUBLE PRECISION :: dt = 0.0, myCommTime, myCompTime,myTime,myTimeMin,myTimeMean
  LOGICAL :: arraySyntax = .false.

  LOGICAL :: uberComm = .false.
  INTEGER :: uberGhost = 15
  
  ! MPI
  CALL MPI_INIT(mpierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

  nargs = command_argument_count()

  ! Print usage:
  IF ( rank == 0 .AND. nargs == 0) THEN
     PRINT*,"USAGE:"
     PRINT*,"Serial: ./miniApp [interations=100,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Parallel: mpirun -n [num procs] miniApp [interations=100,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Examples:"
     PRINT*,"./miniApp 100 32 1 32 1 32 1"
     PRINT*,"mpirun -n 8 ./miniApp 100 64 2 64 2 64 2"
  ENDIF

  ! Default domain, grid and processors map
  x1 = 0.0
  xn = 1.0
  y1 = 0.0
  yn = 1.0
  z1 = 0.0
  zn = 1.0

  ! Grid1
  nx = 32
  ny = 1
  nz = 1

  ! Proc map
  px = 1
  py = 1
  pz = 1

  ! N-species
  ns = 1

  ! Iterations
  iterations = 100

  ! Parse the simple input
  ii=1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') iterations

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') nx

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') px

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') ny

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') py

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') nz

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') pz

  ii=ii+1
  IF (nargs >= ii) CALL GETARG(ii,arg)
  IF (nargs >= ii) READ(arg,'(I10)') ns


  periodic = .false.

  if (.not. periodic) then  
     bx1 = "NONE"
     bxn = "NONE"
     by1 = "NONE"
     byn = "NONE"
     bz1 = "NONE"
     bzn = "NONE"
  else
     bx1 = "PERI"
     bxn = "PERI"
     by1 = "PERI"
     byn = "PERI"
     bz1 = "PERI"
     bzn = "PERI"
  end if

  
  simtime = 0.0D0

  ! Setup matrices/solvers
  CALL setup(0,0,MPI_COMM_WORLD,nx,ny,nz,px,py,pz,0,x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn)
  CALL setup_mesh(0,0)
  CALL point_to_objects(0,0)

  ax = nx / px
  ay = ny / py
  az = nz / pz
  
  ! Allocated some arrays
  ALLOCATE( rho(ax,ay,az) )
  ALLOCATE( u(ax,ay,az) )
  ALLOCATE( v(ax,ay,az) )
  ALLOCATE( w(ax,ay,az) )
  ALLOCATE( et(ax,ay,az) )
  ALLOCATE( ie(ax,ay,az) )
  ALLOCATE( p(ax,ay,az) )
  ALLOCATE( T(ax,ay,az) )
  ALLOCATE( Y(ax,ay,az,ns) )
  

  ALLOCATE( rad(ax,ay,az) )
  ALLOCATE( Fx(ax,ay,az) )
  ALLOCATE( Fy(ax,ay,az) )
  ALLOCATE( Fz(ax,ay,az) )

  ALLOCATE( Fxx(ax,ay,az) )
  ALLOCATE( Fyx(ax,ay,az) )
  ALLOCATE( Fzx(ax,ay,az) )
  ALLOCATE( Fxy(ax,ay,az) )
  ALLOCATE( Fyy(ax,ay,az) )
  ALLOCATE( Fzy(ax,ay,az) )
  ALLOCATE( Fxz(ax,ay,az) )
  ALLOCATE( Fyz(ax,ay,az) )
  ALLOCATE( Fzz(ax,ay,az) )

  
  ALLOCATE( tx(ax,ay,az) )
  ALLOCATE( ty(ax,ay,az) )
  ALLOCATE( tz(ax,ay,az) )
  
  ALLOCATE(uxx(ax,ay,az) )
  ALLOCATE(uyy(ax,ay,az) )
  ALLOCATE(uzz(ax,ay,az) )
  ALLOCATE(uxy(ax,ay,az) )
  ALLOCATE(uxz(ax,ay,az) )
  ALLOCATE(uyz(ax,ay,az) )
  ALLOCATE(uyx(ax,ay,az) )
  ALLOCATE(uzx(ax,ay,az) )
  ALLOCATE(uzy(ax,ay,az) )

  ALLOCATE(beta(ax,ay,az) )
  ALLOCATE(mu(ax,ay,az) )
  ALLOCATE(kappa(ax,ay,az) )
  
  
  ALLOCATE( tmp(ax,ay,az) )
  ALLOCATE( bar(ax,ay,az) )

  
  ALLOCATE( RHS(ax,ay,az,ns+4) )

  if (uberComm) then
     ALLOCATE( RHSp(-uberGhost+1:ax+uberGhost,-uberGhost+1:ay+uberGhost,-uberGhost+1:az+uberGhost, ns+4) )
  end if

  
  ! From is whatebver you want to come back from the device (to host)... impl. does allocate
  ! alloc is only data on the device
  
  
  ! Initialize some profiles
  ! rho = x
  rad = SQRT( ( mesh_ptr%xgrid - 0.5 )**2 + ( mesh_ptr%ygrid - 0.5 )**2 + ( mesh_ptr%zgrid - 0.5 )**2 )
  rho = (1.0 - TANH( (rad - .25 ) / .05 ))*0.5 + 1.0
  u = 0.0
  v = 0.0
  w = 0.0
  ie = 1.0
  Y = 1.0

  CALL EOS(ie,rho,p,t)
  
  ! From(device to host) To(host to device)
  !$omp target data map(to:RHS,u,v,w,Y,p,rho,et,ie,T) &
  !$omp             map(to:mesh_ptr%GridLen) &
  !$omp             map(alloc:Fx,Fy,Fz,Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz,tmp,bar) &
  !$omp             map(alloc:uxx,uyy,uzz,uxy,uxz,uyz,uyx,uzx,uzy) &
  !$omp             map(alloc:beta,mu,kappa,tx,ty,tz)
  
  
  ! Time the derivatives
  CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  call cpu_time(t1c)
  DO tt=1,iterations


     if (uberComm) then
        do eom=1,4+ns
           CALL ghostx(1,RHSp(:,:,:,eom))
           CALL ghosty(1,RHSp(:,:,:,eom))
           CALL ghostz(1,RHSp(:,:,:,eom))
        end do
     end if
        
     
     CALL startCPU()
     !$omp target teams distribute parallel do collapse(3)
     DO j=1,ay
        DO k=1,az
           DO i=1,ax
              ie(i,j,k) = et(i,j,k) - .5 * rho(i,j,k) * &
                   (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     CALL endCPU()
     
     CALL EOS(ie,rho,p,t)
     !CALL EOS_nx(ie,rho,p,t,ax,ay,az)
     
     ! Mass equation
     DO n=1,ns
        CALL startCPU()
        !$omp target teams distribute parallel do collapse(3)
        DO k=1,az           
           DO j=1,ay
              DO i=1,ax
                 Fx(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * u(i,j,k)
                 Fy(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * v(i,j,k)
                 Fz(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * w(i,j,k)
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
        CALL endCPU()

        
        CALL ddx(Fx,tmp,ax,ay,az)
        CALL ddy(Fy,Fx,ax,ay,az)
        !$omp target teams distribute parallel do collapse(3)
        DO k=1,az     
           DO j=1,ay
              DO i=1,ax
                 Fx(i,j,k) = Fx(i,j,k) + tmp(i,j,k)
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
        CALL ddz(Fz,Fy,ax,ay,az)
        !$omp target teams distribute parallel do collapse(3)
        DO k=1,az
           DO j=1,ay
              DO i=1,ax
                 RHS(i,j,k,4+ns) = Fx(i,j,k)+ Fy(i,j,k)
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
     END DO
     
     ! Form the artificial viscosities
     CALL ddx(u,uxx,ax,ay,az)
     CALL ddy(v,uyy,ax,ay,az)
     CALL ddz(w,uzz,ax,ay,az)
     CALL ddy(u,uxy,ax,ay,az)
     CALL ddz(u,uxz,ax,ay,az)
     CALL ddz(v,uyz,ax,ay,az)
     CALL ddx(v,uyx,ax,ay,az)
     CALL ddx(w,uzx,ax,ay,az)
     CALL ddy(w,uzy,ax,ay,az)

     ! Shocks
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = uxx(i,j,k) + uyy(i,j,k) + uzz(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     ! beta = Cb * rho abs( dd4( divU) * del^2 )*
     CALL dd4x( Fx, Fy, ax,ay,az)

     CALL startCPU()
     
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = ABS(Fy(i,j,k)) * mesh_ptr%GridLen(i,j,k)**2
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL endCPU()
     
     CALL gFilter( Fx, Fy, ax,ay,az)

     CALL startCPU()
     
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              beta(i,j,k) = rho(i,j,k) * Fy(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL endCPU()
     
     CALL startCPU()
     

     !$omp target teams distribute parallel do collapse(3) &
     !$omp       private(divu,txx,tyy,tzz,txy,tyz,txz)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax

              divu = uxx(i,j,k) + uyy(i,j,k) + uzz(i,j,k)
              
              txx = p(i,j,k) + beta(i,j,k)*divu + mu(i,j,k)*(uxx(i,j,k)+uxx(i,j,k)-2./3.*divu)
              tyy = p(i,j,k) + beta(i,j,k)*divu + mu(i,j,k)*(uyy(i,j,k)+uyy(i,j,k)-2./3.*divu)
              tzz = p(i,j,k) + beta(i,j,k)*divu + mu(i,j,k)*(uzz(i,j,k)+uzz(i,j,k)-2./3.*divu)

              txy = mu(i,j,k)*(uxy(i,j,k)+uyx(i,j,k) )
              tyz = mu(i,j,k)*(uyz(i,j,k)+uzy(i,j,k) )
              txz = mu(i,j,k)*(uxz(i,j,k)+uzx(i,j,k) )              
              
              
              ! Momentum equation (x)
              Fxx(i,j,k) = rho(i,j,k) * u(i,j,k) * u(i,j,k) + txx
              Fyx(i,j,k) = rho(i,j,k) * u(i,j,k) * v(i,j,k) + txy
              Fzx(i,j,k) = rho(i,j,k) * u(i,j,k) * w(i,j,k) + txz
     
              ! Momentum equation (y)
              Fxy(i,j,k) = rho(i,j,k) * v(i,j,k) * u(i,j,k) + txy
              Fyy(i,j,k) = rho(i,j,k) * v(i,j,k) * v(i,j,k) + tyy
              Fzy(i,j,k) = rho(i,j,k) * v(i,j,k) * w(i,j,k) + tyz
     
              ! Momentum equation (z)
              Fxz(i,j,k) = rho(i,j,k) * w(i,j,k) * u(i,j,k) + txz
              Fyz(i,j,k) = rho(i,j,k) * w(i,j,k) * v(i,j,k) + tyz
              Fzz(i,j,k) = rho(i,j,k) * w(i,j,k) * w(i,j,k) + tzz
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do

     CALL endCPU()
     
     !CALL div(Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz, &
     !RHS(:,:,:,1),RHS(:,:,:,2),RHS(:,:,:,3) )
     
     CALL ddx(Fxx,tmp,ax,ay,az)
     CALL ddy(Fxy,Fx,ax,ay,az)
     !Fx = Fx + Fz
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az     
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = Fx(i,j,k) + tmp(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL ddz(Fxz,Fy,ax,ay,az)
     !RHS(:,:,:,1) = Fx + Fy
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              RHS(i,j,k,1) = Fx(i,j,k)+ Fy(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do

     CALL ddx(Fyx,tmp,ax,ay,az)
     CALL ddy(Fyy,Fx,ax,ay,az)
     !Fx = Fx + Fz
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az     
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = Fx(i,j,k) + tmp(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL ddz(Fyz,Fy,ax,ay,az)
     !RHS(:,:,:,2) = Fx + Fy
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              RHS(i,j,k,2) = Fx(i,j,k)+ Fy(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do

     CALL ddx(Fzx,tmp,ax,ay,az)
     CALL ddy(Fzy,Fx,ax,ay,az)
     !Fx = Fx + Fz
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az     
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = Fx(i,j,k) + tmp(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL ddz(Fzz,Fy,ax,ay,az)
     !RHS(:,:,:,3) = Fx + Fy
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              RHS(i,j,k,3) = Fx(i,j,k)+ Fy(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do

     CALL startCPU()
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              ! Energy equation
              et(i,j,k) = ie(i,j,k) + .5 * rho(i,j,k) * &
                   & (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
              T(i,j,k) = et(i,j,k) * 1.0
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL endCPU()

     
     !CALL grad(T,tx,ty,tz)
     CALL ddx( T, tx, ax, ay, az )
     CALL ddy( T, tx, ax, ay, az )
     CALL ddz( T, tx, ax, ay, az )

     CALL startCPU()
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = et(i,j,k) * u(i,j,k) - tx(i,j,k)
              Fy(i,j,k) = et(i,j,k) * v(i,j,k) - ty(i,j,k)
              Fz(i,j,k) = et(i,j,k) * w(i,j,k) - tz(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     CALL endCPU()
     
     !CALL div(Fx,Fy,Fz,RHS(:,:,:,4))
     CALL ddx(Fx,tmp,ax,ay,az)
     CALL ddy(Fy,Fx,ax,ay,az)
     !Fx = Fx + Fz
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az     
        DO j=1,ay
           DO i=1,ax
              Fx(i,j,k) = Fx(i,j,k) + tmp(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
     
     CALL ddz(Fz,Fy,ax,ay,az)
     !RHS(:,:,:,4) = Fx + Fy
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              RHS(i,j,k,4) = Fx(i,j,k)+ Fy(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do

     ! Integrate the equaions
     CALL startCPU()
     if (arraySyntax) then
        do n=1,ns
           Y(:,:,:,n)  = ( rho*Y(:,:,:,n)  - dt * RHS(:,:,:,4+n)) / rho
        end do
        u  = ( rho*u  - dt * RHS(:,:,:,1))/rho
        v  = ( rho*v  - dt * RHS(:,:,:,2))/rho
        w  = ( rho*w  - dt * RHS(:,:,:,3))/rho
        et = et - dt * RHS(:,:,:,4)
     else

        !$omp target teams distribute parallel do collapse(4) 
        DO n=1,ns
           DO k=1,az
              DO j=1,ay
                 DO i=1,ax                    
                    Y(i,j,k,n) = ( Y(i,j,k,n)*rho(i,j,k) - dt * RHS(i,j,k,4+n))/rho(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
        
        !$omp target teams distribute parallel do collapse(3) private(txx)
        DO k=1,az
           DO j=1,ay
              DO i=1,ax
                 !DO n=1,ns
                 !   Y(i,j,k,n) = ( Y(i,j,k,n)*rho(i,j,k) - dt * RHS(i,j,k,4+n))/rho(i,j,k)
                 !END DO
                 txx = rho(i,j,k)
                 et(i,j,k) = et(i,j,k) - dt * RHS(i,j,k,4)
                 !u(i,j,k)  = ( rho(i,j,k)*u(i,j,k)  - dt * RHS(i,j,k,1))/rho(i,j,k)
                 !v(i,j,k)  = ( rho(i,j,k)*v(i,j,k)  - dt * RHS(i,j,k,2))/rho(i,j,k)
                 !w(i,j,k)  = ( rho(i,j,k)*w(i,j,k)  - dt * RHS(i,j,k,3))/rho(i,j,k)
                 u(i,j,k)  = ( txx*u(i,j,k)  - dt * RHS(i,j,k,1))/txx
                 v(i,j,k)  = ( txx*v(i,j,k)  - dt * RHS(i,j,k,2))/txx
                 w(i,j,k)  = ( txx*w(i,j,k)  - dt * RHS(i,j,k,3))/txx
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
     endif
     CALL endCPU()
     
     ! Filter the equations
     DO n=1,ns
        !$omp target teams distribute parallel do collapse(3)
        DO k=1,az
           DO j=1,ay
              DO i=1,ax
                 tmp(i,j,k) = rho(i,j,k)*Y(i,j,k,n)
              END DO
           END DO
        END DO
        !$omp end target teams distribute parallel do
        
        CALL sFilter( tmp,RHS(:,:,:,n), ax,ay,az)                
        
     END DO

     CALL startCPU()
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              rho(i,j,k) = 0.0D0
              DO n=1,ns
                 rho(i,j,k) = rho(i,j,k) + RHS(i,j,k,n)
              END DO
              DO n=1,ns
                 Y(i,j,k,n) = RHS(i,j,k,n) / rho(i,j,k)
              END DO
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do     
     !DO n=1,ns
     !   Y(:,:,:,n) = rhs(:,:,:,n)/rho
     !END DO

     
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              tmp(i,j,k) = rho(i,j,k)*u(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do
        
     CALL endCPU()
    

     CALL sFilter( tmp, bar, ax,ay,az)

     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              u(i,j,k) = bar(i,j,k) / rho(i,j,k)
              tmp(i,j,k) = rho(i,j,k)*v(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do    
     !u = bar / rho
     !tmp = rho*v

     
     CALL sFilter( tmp, bar, ax,ay,az)
     
     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              v(i,j,k) = bar(i,j,k) / rho(i,j,k)
              tmp(i,j,k) = rho(i,j,k)*w(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do    
     !v = bar / rho     
     !tmp = rho*w


     CALL sFilter( tmp, bar, ax,ay,az)


     !$omp target teams distribute parallel do collapse(3)
     DO k=1,az
        DO j=1,ay
           DO i=1,ax
              w(i,j,k) = bar(i,j,k) / rho(i,j,k)
              tmp(i,j,k) = et(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target teams distribute parallel do 
     !w = bar / rho     
     !tmp = et     
     
     CALL sFilter( tmp, et, ax,ay,az)

     
  END DO

  !$omp end target data

  
  !CALL SYSTEM_CLOCK( t2, clock_rate, clock_max)
  CALL cpu_time( t2c )


  IF ( rank == 0 ) THEN
     print*,'Total time = ', real(t2c-t1c)
  END IF

  
  CALL MPI_ALLREDUCE( comm_time   , myCommTime   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
  CALL MPI_ALLREDUCE( real(t2c-t1c), myCompTime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)

  myCompTime = t2c-t1c - myCommTime
  
  IF ( rank == 0 ) THEN
     print*,'Comm time = ', real(myCommTime) 
  END IF

  IF ( rank == 0 ) THEN
     print*,'CPU time = ', real(myCompTime)
  END IF
  
  IF ( rank == 0 ) THEN
     print*,'COMM/Total time = ', real( myCommTime / (t2c-t1c) )
  END IF

  IF ( rank == 0 ) THEN
     print*,'CPU/Total time = ', real(myCompTime / (t2c-t1c) )
  END IF
  
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, procs)

  DO n=1,6
     CALL MPI_ALLREDUCE( custom_time(n)   , myTime   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
     CALL MPI_ALLREDUCE( custom_time(n)   , myTimeMin   , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
     CALL MPI_ALLREDUCE( custom_time(n)   , myTimeMean   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
    
     IF ( rank == 0 ) THEN
        print*,'Custom = ', real(myTime), real(myTimeMin),real(myTimeMean/procs)
     END IF
  END DO
  
  
  
  CALL remove_objects(0,0)
  CALL MPI_FINALIZE(mpierr)


  CONTAINS
    SUBROUTINE EOS(ie,rho,p,T)
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: ie,rho
      DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: p,t
    END SUBROUTINE EOS
  
    
  
END PROGRAM miniApp




SUBROUTINE EOS(ie,rho,p,T)
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: p,t
  DOUBLE PRECISION :: gamma = 1.4
  INTEGER :: i,j,k

  !$omp target teams distribute parallel do collapse(3) private(gamma)
  DO k=1,size(p,3)
     DO j=1,size(p,2)
        DO i=1,size(p,1)
           p(i,j,k)  = ie(i,j,k)  / rho(i,j,k)  * (gamma - 1.0 )
           t(i,j,k)  = ie(i,j,k)  * (gamma )
        END DO
     END DO
  END DO
  !$omp end target teams distribute parallel do
  
END SUBROUTINE EOS


SUBROUTINE EOS_nx(ie,rho,p,T,nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(OUT) :: p,t
  DOUBLE PRECISION :: gamma = 1.4

  
  p = ie / rho * (gamma - 1.0 )
  t = ie * (gamma )
  
  
END SUBROUTINE EOS_nx
