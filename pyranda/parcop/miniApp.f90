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
  CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn
  REAL(c_double)                 :: simtime
  INTEGER(c_int)                 :: world_id,world_np,mpierr
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz
  REAL(c_double), DIMENSION(:,:,:,:), ALLOCATABLE :: RHS
  INTEGER :: i
  INTEGER :: t1,t2,clock_rate,clock_max
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,ii,iterations
  INTEGER :: rank,ierror
  DOUBLE PRECISION :: dt = 0.0
  !$DEF-FEXL
  
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

  
  ALLOCATE( RHS(ax,ay,az,5) )

  ! Initialize some profiles
  ! rho = x
  rad = SQRT( ( mesh_ptr%xgrid - 0.5 )**2 + ( mesh_ptr%ygrid - 0.5 )**2 + ( mesh_ptr%zgrid - 0.5 )**2 )
  rho = (1.0 - TANH( (rad - .25 ) / .05 ))*0.5 + 1.0
  u = 0.0
  v = 0.0
  w = 0.0
  ie = 1.0

  CALL EOS(ie,rho,p,t)
  
  !$FEXL { vars_to:"rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar,
  !$FEXL            Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz,
  !$FEXL            mesh_ptr%GridLen,RHS" }
  !$END FEXL


  ! Time the derivatives
  CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  DO i=1,iterations
     
     !$FEXL {dim:3,var:['ie','et','rho','u','v','w']}
     ie = et - .5 * rho * (u*u + v*v + w*w )
     !$END FEXL
     
     CALL EOS(ie,rho,p,t)
     
     ! Mass equation

     !$FEXL {dim:3,var:['Fx','Fy','Fz','rho','u','v','w']}
     Fx = rho * u
     Fy = rho * v
     Fz = rho * w
     !$END FEXL     
     CALL div(Fx,Fy,Fz,RHS(:,:,:,1),patch_ptr%ax,patch_ptr%ay,patch_ptr%az)

     ! Momentum equation (x)
     !$FEXL {dim:3,var:['Fxx','Fyx','Fzx','Fxy','Fyy','Fzy','Fxz','Fyz','Fzz','rho','u','v','w','p']}
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
     !$END FEXL
     
     CALL div(Fxx,Fxy,Fxz,Fyx,Fyy,Fyz,Fzx,Fzy,Fzz, &
          RHS(:,:,:,2),RHS(:,:,:,3),RHS(:,:,:,4) )

     ! Energy equation
     !$FEXL {dim:3,var:['et','ie','rho','u','v','w']}
     et = ie + .5 * rho * (u*u + v*v + w*w )
     !$END FEXL
     CALL grad(T,tx,ty,tz)

     !$FEXL {dim:3,var:['Fx','Fy','Fz','et','u','v','w','tx','ty','tz']}
     Fx = et * u - tx
     Fy = et * v - ty
     Fz = et * w - tz
     !$END FEXL
     CALL div(Fx,Fy,Fz,RHS(:,:,:,5),patch_ptr%ax,patch_ptr%ay,patch_ptr%az)

     ! Integrate the equaions
     !$FEXL {dim:3,var:['rho','RHS','u','v','w','et','Fx','Fy','Fz']}
     Fx = rho*u
     Fy = rho*v
     Fz = rho*w
     rho = rho - dt * RHS(:,:,:,1)
     et = et - dt * RHS(:,:,:,5)
     u = (Fx - dt*RHS(:,:,:,2)) / rho
     v = (Fy - dt*RHS(:,:,:,3)) / rho
     w = (Fz - dt*RHS(:,:,:,4)) / rho     
     !$END FEXL
     
     ! Filter the equations
     !$FEXL {dim:3,var:['tmp','rho']}
     tmp = rho
     !$END FEXL
     CALL filter('spectral',tmp,rho)
     !$FEXL {dim:3,var:['tmp','rho','u']}
     tmp = rho*u
     !$END FEXL
     CALL filter('spectral',tmp,bar)
     !$FEXL {dim:3,var:['u','bar','v','tmp','rho']}
     u = bar / rho    
     tmp = rho*v
     !$END FEXL
     CALL filter('spectral',tmp,bar)
     !$FEXL {dim:3,var:['w','bar','v','tmp','rho']}
     v = bar / rho     
     tmp = rho*w
     !$END FEXL
     CALL filter('spectral',tmp,bar)
     !$FEXL {dim:3,var:['w','bar','et','tmp','rho']}
     w = bar / rho     
     tmp = et
     !$END FEXL
     CALL filter('spectral',tmp,et)

     
  END DO


  !$FEXL { vars_from:"rho,u,v,w,et,p,rad,T,ie,Fx,Fy,Fz,tx,ty,tz,tmp,bar,
  !$FEXL              Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz,
  !$FEXL              mesh_ptr%GridLen,RHS" }
  !$END FEXL
  
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
  !$DEF-FEXL
  
  !$FEXL {dim:3,var:['p','ie','rho','t']}
  p = ie / rho * (gamma - 1.0 )
  t = ie * (gamma )
  !$END FEXL
  
END SUBROUTINE EOS


SUBROUTINE EOS_nx(ie,rho,p,T,nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(OUT) :: p,t
  DOUBLE PRECISION :: gamma = 1.4

  
  p = ie / rho * (gamma - 1.0 )
  t = ie * (gamma )
  
  
END SUBROUTINE EOS_nx
