PROGRAM test_laplace

  USE iso_c_binding
  USE LES_compact, ONLY : compact_type
  USE LES_patch, ONLY : patch_type
  USE LES_comm, ONLY : comm_type, LES_comm_world
  USE LES_mesh, ONLY : mesh_type
  USE LES_objects 
  USE parcop, ONLY : setup,plaplacian
  USE nvtx
  !USE cudafor
  IMPLICIT NONE
  INCLUDE "mpif.h"

  
  INTEGER(c_int)                 :: nx,ny,nz,px,py,pz
  REAL(c_double)                 :: x1,xn,y1,yn,z1,zn
  CHARACTER(KIND=c_char,LEN=4)   :: bx1,bxn,by1,byn,bz1,bzn
  REAL(c_double)                 :: simtime
  INTEGER(c_int)                 :: world_id,world_np,mpierr
  REAL(c_double), DIMENSION(:,:,:), ALLOCATABLE :: rho,drho
  INTEGER :: i 
  INTEGER :: t1,t2,clock_rate,clock_max
  CHARACTER(LEN=32) :: arg
  INTEGER :: nargs,ii,iterations
  INTEGER :: rank,ierror

  ! MPI
  CALL MPI_INIT(mpierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

  nargs = command_argument_count()
  
  ! Print usage:
  IF ( rank == 0 .AND. nargs == 0) THEN
     PRINT*,"USAGE:"
     PRINT*,"Serial: ./test_laplace [interations=100,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Parallel: mpirun -n [num procs] test_laplace [interations=100,nx=32,px=1,ny=1,py=1,nz=1,pz=1]"
     PRINT*,"Examples:"
     PRINT*,"./test_laplace 100 32 1 32 1 32 1"
     PRINT*,"mpirun -n 8 ./test_laplace 100 64 2 64 2 64 2"
  ENDIF
  
  ! Default domain, grid and processors map
  x1 = 0.0
  xn = 1.0
  y1 = 0.0
  yn = 1.0
  z1 = 0.0
  zn = 1.0

  ! Grid
  nx = 32
  ny = 1
  nz = 1

  ! Proc map
  px = 1
  py = 1
  pz = 1
  
  ! Iterations
  iterations = 4

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
  CALL setup(0,0,MPI_COMM_WORLD,nx,ny,nz,px,py,pz,x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn)

  !i = cudaMemAdvise(rho,nx/px*ny/py*nz/pz, cudaMemAdviseSetPreferredLocation, 0)
  !i = cudaMemAdvise(drho,nx/px*ny/py*nz/pz, cudaMemAdviseSetPreferredLocation, 0)
  ! Allocated some arrays
  ALLOCATE( rho(nx/px,ny/py,nz/pz) )
  ALLOCATE(drho(nx/px,ny/py,nz/pz) )
  
  ! rho = x
  rho = mesh_ptr%xgrid

  ! Time the derivatives
  CALL SYSTEM_CLOCK( t1, clock_rate, clock_max)
  !$acc data copy(rho) copyout(drho)
  DO i=1,iterations
     PRINT *, "iter ", i, " of ", iterations
     CALL nvtxStartRange("plaplacian")
     CALL plaplacian(rho,drho,nx/px,ny/py,nz/pz)    
     CALL nvtxEndRange()
     IF (.true.) THEN
        !$acc kernels
        rho(:,:,:) = rho(:,:,:)-0.0001*drho(:,:,:)
        !$acc end kernels
        PRINT *, "rho(31, 31, 21) = ", rho(31,31,21)
     ENDIF
  END DO
  !$acc end data
  CALL SYSTEM_CLOCK( t2, clock_rate, clock_max)


  IF ( rank == 0 ) THEN
     print*,'Elapsed time = ', real(t2-t1) / real(clock_rate)
  END IF
  
  CALL remove_objects(0,0)

END PROGRAM test_laplace


