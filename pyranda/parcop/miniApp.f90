PROGRAM miniApp

  USE iso_c_binding
  USE LES_compact, ONLY : compact_type
  USE LES_patch, ONLY : patch_type
  USE LES_comm, ONLY : comm_type, LES_comm_world
  USE LES_mesh, ONLY : mesh_type
  USE LES_objects
  USE parcop, ONLY : setup,ddx,point_to_objects,setup_mesh,grad,filter,div
#ifdef fexlpool
  USE LES_ompsync
#endif
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
  REAL(c_double) :: t1,t2,clock_rate,clock_max,t3,t4
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,ii,iterations
  INTEGER :: rank,ierror
  DOUBLE PRECISION :: dt = 0.0, memtime = 0.0
  DOUBLE PRECISION :: answer = 830909.7500
  integer :: ifunr,jfunr,kfunr,nfunr 

  
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
  u = 1.0
  v = 0.0
  w = 0.0
  ie = 1.0

  CALL EOS(ie,rho,p,t,patch_ptr%ax,patch_ptr%ay,patch_ptr%az)
  
  ! Persistent Data Map.. should always be on the device
  !$omp target data map(tofrom:rho,u,v,w,et,p,rad,T,ie) &
  !$omp             map(tofrom:mesh_ptr%GridLen,RHS)

#ifdef fexlpool
  !! FEXL-POOL version
  CALL fexlPool_setPoolDepth(25,ax,ay,az)
#endif
  
  ! Time the derivatives
  !CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  CALL CPU_TIME(t1)
  
  DO i=1,iterations


     !CALL SYSTEM_CLOCK( t3, clock_rate, clock_max )
     CALL CPU_TIME(t3)
     
#ifdef omppool
     ! Temp. memory
     !$omp target data map(alloc:Fx,Fy,Fz,tx,ty,tz,tmp,bar,Fxx,Fyx,Fzx,Fxy,Fyy,Fzy,Fxz,Fyz,Fzz)
#endif
     
#ifdef fexlpool
     ASSOCIATE(Fx=>fexlPool_get(ax,ay,az),Fy=>fexlPool_get(ax,ay,az),Fz=>fexlPool_get(ax,ay,az) &
          ,tx=>fexlPool_get(ax,ay,az),ty=>fexlPool_get(ax,ay,az),tz=>fexlPool_get(ax,ay,az) &
          ,tmp=>fexlPool_get(ax,ay,az),bar=>fexlPool_get(ax,ay,az),Fxx=>fexlPool_get(ax,ay,az) &
          ,Fyx=>fexlPool_get(ax,ay,az),Fzx=>fexlPool_get(ax,ay,az),Fxy=>fexlPool_get(ax,ay,az) &
          ,Fyy=>fexlPool_get(ax,ay,az),Fzy=>fexlPool_get(ax,ay,az),Fxz=>fexlPool_get(ax,ay,az) &
          ,Fyz=>fexlPool_get(ax,ay,az),Fzz=>fexlPool_get(ax,ay,az)  )
#endif
     

       !CALL SYSTEM_CLOCK( t4, clock_rate, clock_max )
     CALL CPU_TIME( t4 )
       

     !memtime = memtime + real(t4-t3) / real(clock_rate)
     memtime = memtime + (t4-t3)
     

!$omp target teams distribute parallel do collapse(3)
     do kfunr=1,size(ie,3)
       do jfunr=1,size(ie,2)
         do ifunr=1,size(ie,1)
           ie(ifunr,jfunr,kfunr)= et(ifunr,jfunr,kfunr)- .5*rho(ifunr,jfunr,kfunr)*(u(ifunr,jfunr,kfunr)* u(ifunr,jfunr,kfunr)+ v(ifunr,jfunr,kfunr)* v(ifunr,jfunr,kfunr)+ w(ifunr,jfunr,kfunr)* w(ifunr,jfunr,kfunr))
         end do
       end do
     end do
!$omp end target teams distribute parallel do

     
     CALL EOS(ie,rho,p,t,patch_ptr%ax,patch_ptr%ay,patch_ptr%az)

     ! Rest u

!$omp target teams distribute parallel do collapse(3)
     do kfunr=1,size(u,3)
       do jfunr=1,size(u,2)
         do ifunr=1,size(u,1)
           u(ifunr,jfunr,kfunr)= 1.0D0
         end do
       end do
     end do
!$omp end target teams distribute parallel do

     
     ! Mass equation


!$omp target teams distribute parallel do collapse(3)
     do kfunr=1,size(Fx,3)
       do jfunr=1,size(Fx,2)
         do ifunr=1,size(Fx,1)
           Fx(ifunr,jfunr,kfunr)= rho(ifunr,jfunr,kfunr)* u(ifunr,jfunr,kfunr)
            Fy(ifunr,jfunr,kfunr)= rho(ifunr,jfunr,kfunr)* v(ifunr,jfunr,kfunr)
            Fz(ifunr,jfunr,kfunr)= rho(ifunr,jfunr,kfunr)* w(ifunr,jfunr,kfunr)
          end do
       end do
     end do
!$omp end target teams distribute parallel do

     CALL div(Fx,Fy,Fz,RHS(:,:,:,1),patch_ptr%ax,patch_ptr%ay,patch_ptr%az)


!$omp target teams distribute parallel do collapse(3)
     do kfunr=1,size(u,3)
       do jfunr=1,size(u,2)
         do ifunr=1,size(u,1)
           u(ifunr,jfunr,kfunr)= u(ifunr,jfunr,kfunr)* 0.0+RHS(ifunr,jfunr,kfunr, 1 )
         end do
       end do
     end do
!$omp end target teams distribute parallel do

     
#ifdef fexlpool
   END ASSOCIATE
   ierror = fexlPool_free(17)
#endif

#ifdef omppool
   !$omp end target data 
#endif
   
  END DO

  !CALL SYSTEM_CLOCK( t2, clock_rate, clock_max)
  CALL CPU_TIME( t2 )
  
  !$omp end target data  

  IF ( rank == 0 ) THEN
     print*,'Ellapsed time = ', real(t2-t1) !/ real(clock_rate)
     print*,'Memory time = ', real(memtime)
     print*,'Result = ' , real( SUM(ABS(u)) )
     print*,'Answer = ' , real( answer )
  END IF

  CALL remove_objects(0,0)
  CALL MPI_FINALIZE(mpierr)


    
  
END PROGRAM miniApp




SUBROUTINE EOS(ie,rho,p,T,nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(IN) :: ie,rho
  DOUBLE PRECISION, DIMENSION(nx,ny,nz), INTENT(OUT) :: p,t
  DOUBLE PRECISION, DIMENSION(nx,ny,nz) :: tmp1,tmp2
  DOUBLE PRECISION :: gamma = 1.4

  !$omp target data map(alloc:tmp1,tmp2)


!$omp target teams distribute parallel do collapse(3)
  do kfunr=1,size(p,3)
    do jfunr=1,size(p,2)
      do ifunr=1,size(p,1)
        p(ifunr,jfunr,kfunr)= ie(ifunr,jfunr,kfunr)/ rho(ifunr,jfunr,kfunr)*(gamma-1.0 )
        t(ifunr,jfunr,kfunr)= ie(ifunr,jfunr,kfunr)*(gamma )
      end do
    end do
  end do
!$omp end target teams distribute parallel do

  
  !$omp end target data


END SUBROUTINE EOS
