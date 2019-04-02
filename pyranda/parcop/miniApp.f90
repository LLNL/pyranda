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

  
  INTEGER(c_int)                 :: nx,ny,nz,px,py,pz,ax,ay,az
  REAL(c_double)                 :: x1,xn,y1,yn,z1,zn
  INTEGER(c_int)                 :: ns
  CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn
  REAL(c_double)                 :: simtime
  INTEGER(c_int)                 :: world_id,world_np,mpierr
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz
  REAL(c_double), DIMENSION(:,:,:,:), ALLOCATABLE :: RHS,Y
  INTEGER :: i,n
  INTEGER :: t1,t2,clock_rate,clock_max
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,io,iterations
  INTEGER :: rank,ierror
  DOUBLE PRECISION :: dt = 0.0
  !$DEF
  
  ! MPI
  CALL MPI_INIT(mpierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

  nargs = command_argument_count()

  ! Print usage:
  IF ( rank == 0 .AND. nargs == 0) THEN
     PRINT*,"USAGE:"
     PRINT*,"Serial: ./miniApp [interations=100,ns=3,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Parallel: mpirun -n [num procs] miniApp [interations=100,ns=3,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Examples:"
     PRINT*,"./miniApp 100 3 32 1 32 1 32 1"
     PRINT*,"mpirun -n 8 ./miniApp 100 3 64 2 64 2 64 2"
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

  ! Iterations
  iterations = 100

  ! Parse the simple input
  io=1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') iterations

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') ns

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') nx

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') px

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') ny

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') py

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') nz

  io=io+1
  IF (nargs >= io) CALL GETARG(io,arg)
  IF (nargs >= io) READ(arg,'(I10)') pz



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

  ! NS species array
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

  
  ALLOCATE( RHS(ax,ay,az,4+ns) )

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
  
 

  ! Start the pseudo-physics loop
  CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  
  DO i=1,iterations

     ie = et - .5 * rho * (u*u + v*v + w*w )
     CALL EOS(ie,rho,p,t)
     !CALL EOS_nx(ie,rho,p,t,ax,ay,az)
     
     CALL grad(T,tx,ty,tz)

     
     ! Mass equation(s)
     DO n=1,ns
        Fx = rho * u * Y(:,:,:,n)
        Fy = rho * v * Y(:,:,:,n)
        Fz = rho * w * Y(:,:,:,n)
        CALL div(Fx,Fy,Fz,RHS(:,:,:,4+n))
     END DO


     !$UNROLL {dim:3,var:['Fxx','Fyx','Fzx','rho','u','v','w','p','Fxy','Fyy','Fzy','Fxz','Fyz','Fzz','ie','et','Fx','Fy','Fz','tx','ty','tz']}
     ! Momentum equation (x)
     Fxx = rho * u * u + p
     Fyx = rho * u * v
     Fzx = rho * u * w 
     
     ! Momentum equation (y)
     Fxy = rho * v * u 
     Fyy = rho * v * v + p
     Fzy = rho * v * w 
     
     ! Momentum equation (z)
     Fxz = rho * w * u 
     Fyz = rho * w * v
     Fzz = rho * w * w + p

     ! Energy equation
     et = ie + .5 * rho * (u*u + v*v + w*w)
     
     Fx = et * u - tx
     Fy = et * v - ty
     Fz = et * w - tz
     !$END UNROLL

     CALL div(Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz, &
          RHS(:,:,:,1),RHS(:,:,:,2),RHS(:,:,:,3) )          
     CALL div(Fx,Fy,Fz,RHS(:,:,:,4))

     ! Integrate the equaions
     rho = rho - dt * RHS(:,:,:,1)

     tmp = rho*u - dt * RHS(:,:,:,2)
     u = tmp / rho

     tmp = rho*v - dt * RHS(:,:,:,3)
     v = tmp / rho

     tmp = rho*w - dt * RHS(:,:,:,4)
     w = tmp / rho
          
     et = et - dt * RHS(:,:,:,5)

     
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

  
  p = ie * (gamma - 1.0 )
  t = p / rho 
  
  
END SUBROUTINE EOS


SUBROUTINE EOS_nx(ie,rho,p,T,nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(OUT) :: p,t
  DOUBLE PRECISION :: gamma = 1.4

  
  p = ie * (gamma - 1.0 )
  t = p / rho
  
  
END SUBROUTINE EOS_nx
