PROGRAM miniApp

  USE iso_c_binding
  USE LES_compact, ONLY : compact_type
  USE LES_patch, ONLY : patch_type
  USE LES_comm, ONLY : comm_type, LES_comm_world
  USE LES_mesh, ONLY : mesh_type
  USE LES_objects
  USE parcop, ONLY : setup,ddx,point_to_objects,setup_mesh,grad,filter,div
  IMPLICIT NONE
  INCLUDE "mpif.h"

  
  INTEGER(c_int)                 :: nx,ny,nz,px,py,pz,ax,ay,az,ns
  REAL(c_double)                 :: x1,xn,y1,yn,z1,zn
  CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn
  REAL(c_double)                 :: simtime
  INTEGER(c_int)                 :: world_id,world_np,mpierr
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz
  REAL(c_double), DIMENSION(:,:,:,:), ALLOCATABLE :: RHS,Y
  INTEGER :: tt,i,j,k,n
  INTEGER :: t1,t2,clock_rate,clock_max
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,ii,iterations
  INTEGER :: rank,ierror
  DOUBLE PRECISION :: dt = 0.0
  LOGICAL :: arraySyntax = .false.
  
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



  bx1 = "NONE"
  bxn = "NONE"
  by1 = "NONE"
  byn = "NONE"
  bz1 = "NONE"
  bzn = "NONE"

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

  ALLOCATE( tmp(ax,ay,az) )
  ALLOCATE( bar(ax,ay,az) )

  
  ALLOCATE( RHS(ax,ay,az,ns+4) )

  ! From is whatebver you want to come back from the device (to host)... impl. does allocate
  ! alloc is only data on the device
  
  !$omp target data map(from:RHS,u,v,w,Y,p,rho) map(alloc:Fx,Fy,Fz,Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz)
  
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
 

  ! Time the derivatives
  CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  DO tt=1,iterations

     DO i=1,ax
        DO j=1,ay
           DO k=1,az
              ie(i,j,k) = et(i,j,k) - .5 * rho(i,j,k) * &
                   (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
           END DO
        END DO
     END DO
 
     CALL EOS(ie,rho,p,t)
     !CALL EOS_nx(ie,rho,p,t,ax,ay,az)
     
     ! Mass equation
     DO n=1,ns
        !$omp target parallel collapse(3)
        DO i=1,ax
           DO j=1,ay
              DO k=1,az
                 Fx(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * u(i,j,k)
                 Fy(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * v(i,j,k)
                 Fz(i,j,k) = Y(i,j,k,ns) * rho(i,j,k) * w(i,j,k)
              END DO
           END DO
        END DO
        !$end omp target parallel
        CALL div(Fx,Fy,Fz,RHS(:,:,:,4+ns))        
     END DO

     !$omp target parallel collapse(3)
     DO i=1,ax
        DO j=1,ay
           DO k=1,az
              ! Momentum equation (x)
              Fxx(i,j,k) = rho(i,j,k) * u(i,j,k) * u(i,j,k) + p(i,j,k)
              Fyx(i,j,k) = rho(i,j,k) * u(i,j,k) * v(i,j,k)
              Fzx(i,j,k) = rho(i,j,k) * u(i,j,k) * w(i,j,k) 
     
              ! Momentum equation (y)
              Fxy(i,j,k) = rho(i,j,k) * v(i,j,k) * u(i,j,k) 
              Fyy(i,j,k) = rho(i,j,k) * v(i,j,k) * v(i,j,k) + p(i,j,k)
              Fzy(i,j,k) = rho(i,j,k) * v(i,j,k) * w(i,j,k) 
     
              ! Momentum equation (z)
              Fxz(i,j,k) = rho(i,j,k) * w(i,j,k) * u(i,j,k) 
              Fyz(i,j,k) = rho(i,j,k) * w(i,j,k) * v(i,j,k)
              Fzz(i,j,k) = rho(i,j,k) * w(i,j,k) * w(i,j,k) + p(i,j,k)
           END DO
        END DO
     END DO
     !$omp end target parallel
     
     CALL div(Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz, &
          RHS(:,:,:,1),RHS(:,:,:,2),RHS(:,:,:,3) )

     !$omp target parallel collapse(3)
     DO i=1,ax
        DO j=1,ay
           DO k=1,az
              ! Energy equation
              et(i,j,k) = ie(i,j,k) + .5 * rho(i,j,k) * &
                   & (u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k) )
              T(i,j,k) = et(i,j,k) * 1.0
           END DO
        END DO
     END DO
     !$end omp target parallel
     
     CALL grad(T,tx,ty,tz)

     !$omp target parallel collapse(3)
     DO i=1,ax
        DO j=1,ay
           DO k=1,az
              Fx(i,j,k) = et(i,j,k) * u(i,j,k) - tx(i,j,k)
              Fy(i,j,k) = et(i,j,k) * v(i,j,k) - ty(i,j,k)
              Fz(i,j,k) = et(i,j,k) * w(i,j,k) - tz(i,j,k)
           END DO
        END DO
     END DO
     !$end omp target parallel
     
     CALL div(Fx,Fy,Fz,RHS(:,:,:,4))
     
     ! Integrate the equaions
     if (arraySyntax) then
        do n=1,ns
           Y(:,:,:,n)  = ( rho*Y(:,:,:,n)  - dt * RHS(:,:,:,4+n)) / rho
        end do
        u  = ( rho*u  - dt * RHS(:,:,:,1))/rho
        v  = ( rho*v  - dt * RHS(:,:,:,2))/rho
        w  = ( rho*w  - dt * RHS(:,:,:,3))/rho
        et = et - dt * RHS(:,:,:,4)
     else
        !$omp target parallel collapse(3)
        DO i=1,ax
           DO j=1,ay
              DO k=1,az
                 DO n=1,ns
                    Y(i,j,k,ns) = ( Y(i,j,k,n)*rho(i,j,k) - dt * RHS(i,j,k,1))/rho(i,j,k)
                 END DO
                 et(i,j,k) = et(i,j,k) - dt * RHS(i,j,k,5)
                 u(i,j,k)  = ( rho(i,j,k)*u(i,j,k)  - dt * RHS(i,j,k,2))/rho(i,j,k)
                 v(i,j,k)  = ( rho(i,j,k)*v(i,j,k)  - dt * RHS(i,j,k,2))/rho(i,j,k)
                 w(i,j,k)  = ( rho(i,j,k)*w(i,j,k)  - dt * RHS(i,j,k,2))/rho(i,j,k)
              END DO
           END DO
        END DO
        !$end omp target parallel
     endif
     
     ! Filter the equations
     tmp = rho
     CALL filter('spectral',tmp,rho)
     
     tmp = rho*u
     CALL filter('spectral',tmp,bar)
     u = bar / rho
     
     tmp = rho*v
     CALL filter('spectral',tmp,bar)
     v = bar / rho
     
     tmp = rho*w
     CALL filter('spectral',tmp,bar)
     w = bar / rho
     
     tmp = et
     CALL filter('spectral',tmp,et)

     
  END DO

  !$omp end target data

  
  CALL SYSTEM_CLOCK( t2, clock_rate, clock_max)


  IF ( rank == 0 ) THEN
     print*,'Ellapsed time = ', real(t2-t1) / real(clock_rate)
  END IF

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

  !$omp target parallel collapse(3)
  DO i=1,size(p,1)
     DO j=1,size(p,2)
        DO k=1,size(p,3)
           p(i,j,k)  = ie(i,j,k)  / rho(i,j,k)  * (gamma - 1.0 )
           t(i,j,k)  = ie(i,j,k)  * (gamma )
        END DO
     END DO
  END DO
  !$end omp
  
END SUBROUTINE EOS


SUBROUTINE EOS_nx(ie,rho,p,T,nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(OUT) :: p,t
  DOUBLE PRECISION :: gamma = 1.4

  
  p = ie / rho * (gamma - 1.0 )
  t = ie * (gamma )
  
  
END SUBROUTINE EOS_nx
