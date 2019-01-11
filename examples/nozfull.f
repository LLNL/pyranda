! you must re-define this module name to be unique across all problems
MODULE nozfull_data
  USE mpi
  SAVE
  ! ---------------------------------------------------------------------------
  ! Place Problem dependent state data here, similar to inputs.f and globals.f
  ! ---------------------------------------------------------------------------  
  CHARACTER(LEN=30) 	:: init_type    = '2D'
  CHARACTER(LEN=90) 	:: restart_dir  = '/path/to/directory'
  CHARACTER(LEN=90) 	:: gridpath    = 'nozfull.grid'
  CHARACTER(LEN=90) 	:: gridfile    = 'nozfull.grid'
  CHARACTER(LEN=90) 	:: flowfile    = 'nozfull.inflow'
  CHARACTER(LEN=90) 	:: initfile    = 'nozfull.init'
  DOUBLE PRECISION 	:: Pinitial     = 1.01d6 
  DOUBLE PRECISION 	:: Tinitial     = 300.0d0
  DOUBLE PRECISION 	:: Pnew        = 1.0d6
  DOUBLE PRECISION 	:: Pold        = 2.0d6
  LOGICAL               :: Ptrans      = .FALSE.
  DOUBLE PRECISION 	:: Tnew        = 300.0d0
  DOUBLE PRECISION 	:: Told        = 300.0d0
  LOGICAL               :: Ttrans      = .FALSE.
  ! Not set in Namelist, these variables are read in and used to set the inflow BC
  CHARACTER(LEN=30) 	:: prob_jobname = 'nozfull'
  DOUBLE PRECISION 	:: Mach  = 1.0d0 
  DOUBLE PRECISION      :: P_in   = 1.01d6
  DOUBLE PRECISION      :: rho_in = 1.0d0
  DOUBLE PRECISION      :: U_in = 1000.0d0
  DOUBLE PRECISION      :: T_in = 300.0d0
  DOUBLE PRECISION      :: e_in = 1.0d0      
  DOUBLE PRECISION      :: NPR   = 1.5d0
  DOUBLE PRECISION      :: rho_amb   = 1.0D0
  DOUBLE PRECISION      :: e_amb = 1.0d0
  DOUBLE PRECISION      :: Rgas = 1.0D0
  REAL,ALLOCATABLE, DIMENSION(:,:)  :: randBS
  !! For the Recycling rescaling routines
  DOUBLE PRECISION 	:: del_BL       = 1.0d-1 
  DOUBLE PRECISION 	:: del_star     = 1.0d-2 
  DOUBLE PRECISION 	:: Re_BL        = 1.0d4
  DOUBLE PRECISION 	:: BLalpha      = 1.5D0      
  DOUBLE PRECISION 	:: mu_0         = 1.0D0      
  DOUBLE PRECISION 	:: Pr           = 0.7D0    
  INTEGER        	:: rcy_pt_g     = 100  
  DOUBLE PRECISION 	:: Len_BL       = 4.0d0  !! Only used to get BL_alpha based on Urban and Knight... better to set explicitly
  INTEGER               :: nvar         = 4
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Qave
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: y_r,y_G
  DOUBLE PRECISION :: simtime_old,tauX
  INTEGER :: res_dump
  INTEGER :: rand_seed = 12
  LOGICAL :: init_stats = .FALSE.     ! Initialize the temporal averages? Will do this if simtime < dt (at startup)
  LOGICAL :: BL_flag = .TRUE.              ! A flag to call routine once in BC
  CHARACTER(len=80) :: uaveDir

  !! Digital Filter routines
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A11,A22,A33,A12,Umean,Tmean,RHO_u,RHO_v,RHO_w
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: buI,bvI,bwI,buO,bvO,bwO
  INTEGER :: UIx,UIz,VIx,VIz,WIx,WIz
  INTEGER, DIMENSION(2) :: UIy,VIy,WIy    
  INTEGER :: Nbuff

  !! Stuff to save for RR communication
  INTEGER, DIMENSION(:,:,:,:,:), ALLOCATABLE :: MAP
  INTEGER, DIMENSION(:), ALLOCATABLE :: pmapID
  INTEGER :: rstat(MPI_STATUS_SIZE)


  LOGICAL :: pdonor,precv
  INTEGER :: idonor,irecv
  INTEGER :: yproc
  LOGICAL :: flippy




  DOUBLE PRECISION :: tflip,tflow,tiid
  INTEGER :: step_old

  INTEGER :: check = 0
  INTEGER :: filmax = 30
  LOGICAL :: Ifil = .FALSE.

END MODULE nozfull_data

! -----------------------------------------------------------------------------
! prob_inputs
! Read problem-dependent inputs and save in module variables
! must also define 'jobname'
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_inputs(fileName)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE nozfull_data
  USE globals, ONLY : inputUnit
  USE inputs, ONLY : nx,ny,nz
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fileName
  INTEGER :: nxp,nyp
  INTEGER :: in_unit 
  INTEGER, DIMENSION(1) :: seed

  ! Uncomment to define a namelist
  NAMELIST /nozfull_vars/ Tinitial,Pinitial,del_BL,Re_BL,Pr,rcy_pt_g,init_type,restart_dir,init_stats,BL_flag,&
      & gridpath,gridfile,flowfile,initfile, Ifil,filmax,BLalpha

  ! Uncomment to open and read a namelist file
  OPEN(UNIT=inputUnit,FILE=TRIM(fileName),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=inputUnit,NML=nozfull_vars)
  CLOSE(inputUnit)
  
  ! give the right name
  jobname = TRIM(prob_jobname)
  
  ! Read the inflow boundary conditions                                     
  !in_unit = 26
  !OPEN(UNIT=in_unit, file=TRIM(flowfile))
  !READ(in_unit,*) NPR,P_in,rho_in,Mach
  !CLOSE(in_unit)
  NPR = 1.7
  P_in = 0.5283
  rho_in = 0.63395
  Mach = .99997

  ! Set some time independent random numbers
  !ALLOCATE(randBS(2,nz))
  !seed = 12
  !CALL random_seed(put=seed)
  !CALL random_number(randBS)
  
END SUBROUTINE prob_inputs



! -------------------------------------------------------------------------------------
! Set up problem geometry. (This is called on startup and restart before timestepping.)
! -------------------------------------------------------------------------------------
 SUBROUTINE prob_setup()
  USE mpi
  USE constants
  USE inputs
  USE globals
  USE extend
  USE nozfull_data
  USE metrics, ONLY: x_c
  USE interfaces, ONLY : spherecoord, ERF, ghost, SUMproc, ran1,set_time_seed
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Gplane,Lplane
  DOUBLE PRECISION :: delX
  DOUBLE PRECISION, DIMENSION(1) :: randT
  INTEGER :: i,funit
  CHARACTER(LEN=90) :: flipfile,fform
  LOGICAL :: fexist


  ! Some prob constants which are input dependent
  Rgas = Runiv/molwts(1)                        ! Specific Gas constant
  rho_in = rho_in*NPR*Pinitial/(Tinitial*Rgas)  ! Density of inflow- Convert to physical units 
  P_in = P_in*NPR*Pinitial                      ! Pressure of inflow- Convert to physical units
  U_in = Mach*sqrt(P_in*gamma/rho_in)           ! Get the velocity on the inflow
  e_in = (P_in/(gamma-one))/rho_in              ! Energy of inflow- Convert to physical units 
  T_in = P_in / (rho_in * Rgas)                 ! Temp. of inflow- Convert to physical units 
  mu_0 = U_in * rho_in * del_BL / Re_BL         ! Physical viscosity based on inlet parameters
  
  rho_amb = Pinitial/(Tinitial*Rgas)
  e_amb = (Pinitial/(gamma-one))/rho_amb 

  ! Set-up the grid and metric terms for the calculation
  !CALL prob_geom()
  CALL prob_geom_2d()

  !IF (xyzcom_id == master) print*,'test'
  check = 0


  CALL setup_DFinflow()


  !! Get Boundary Layer length from lower right hand corner
  !delX = 0.0D0
  !IF(ix(1)==1 .AND. iy(1)==1 .AND. iz(1)==1) delX = ABS( x_c(2,1,1)-x_c(1,1,1) )
  !delX = SUMproc(delX)
  !Len_BL = dble(rcy_pt_g) * delX


     
END SUBROUTINE prob_setup

! -----------------------------------------------------------------------------
! prob_init
! Initialize state variables at start of run
! Called by the Miranda main program.
! -----------------------------------------------------------------------------
SUBROUTINE prob_init(rho,u,v,w,e,Y,p,T)
  USE mpi
  USE constants
  USE globals, ONLY: simtime,molwts,ax,ay,az,ix,iy,iz,xloc,yloc,zloc,x1proc,xnproc,y1proc
  USE globals, ONLY: iodata,nres,iodir,jobdir
  USE interfaces, ONLY: restart,filter,boundary
  USE inputs
  USE prob_interface, ONLY : jobname
  USE nozfull_data
  USE metrics, ONLY : x_c,y_c
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: rho,u,v,w,e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(OUT) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT)   :: p,T
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dum,tmp0
  DOUBLE PRECISION, DIMENSION(nx,ny) :: tmpP,tmpRho,tmpU,tmpV
  LOGICAL :: noz
  INTEGER :: in_unit,i,j,k,re_unit  
  DOUBLE PRECISION :: Mflow,Pflow,Gam,mwt,thick,x0,y0,jump,dumy
  INTEGER :: nxp,nyp
  CHARACTER(LEN=90) :: comments,filename


      p = Pinitial
      T = Tinitial
      rho = P / (T*Rgas)

      e = (p/(gamma-one))/rho
      u = zero
      v = zero
      w = zero
      Y = one
      


      DO i=1,ax
         IF (ix(i) < 350) THEN
            u(i,:,:) = Umean
            T(i,:,:) = Tmean
            P(i,:,:) = P_in
            rho(i,:,:) = P_in / (Tmean*Rgas) 
            e(i,:,:) = (P_in/(gamma-one))/rho(i,:,:)
         END IF
      END DO

      CALL filter(gfilter,u,u)

    SELECT CASE(init_type)
         CASE('3D')
         !  Another option here to read in a restart file as the initialization of the flow field.
         !  The global variable 'iodata' is a pointer to the global flow field variables.  We only need
         !  to set iodata for a full initialization.  A new restart file will be written with the identical 
         !  data.
         WRITE(iodir,'(2A)') TRIM(jobdir) , '/res_init'
         CALL restart('r',TRIM(iodir),iodata(:,:,:,1:nres))

         DO i=1,nres
            CALL filter(gfilter,iodata(:,:,:,i),iodata(:,:,:,i))
         END DO

     END SELECT



      ! This inits these variables.  Setup will try to read in RHOs from restart dir
      RHO_u = zero
      RHO_v = zero
      RHO_w = zero
      CALL DFinflow(rho,u,v,w,e) 
      

END SUBROUTINE prob_init

! -----------------------------------------------------------------------------
! prob_eos
! Problem specific equation of state call
! Called by routine "eos()" in eos.f in the case when eostype == 'CASES'
! -----------------------------------------------------------------------------
 SUBROUTINE  prob_eos(rho,e,Y,p,T,c,dTdE,kR,kP,dkPdE)
  USE mpi
  USE constants, only : zero
  USE prob_interface, ONLY : jobname
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: rho,e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: p,T
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: c,dTdE,kR,kP,dkPdE
  ! should never be called for Incompressible problem
  CALL fatal(-1,"prob_eos called for Incompressible problem ",jobname)
  ! give INTENT(OUT) args a value to get rid of annoying compiler warnings
  if (present(c)) c = zero
  if (present(dTde)) dTde = zero
  if (present(kR)) kR = zero
  if (present(kP)) kP = zero
  if (present(dkPdE)) dkPdE = zero
 END SUBROUTINE prob_eos
 
! ----------------------------------------------------------------------------------
! prob_stats
! Problem specific function to compute and write statistics data every "stats"
! timesteps
! Called by routine "monitor()".
! -----------------------------------------------------------------------------------
SUBROUTINE prob_stats(simtime,rho,u,v,w,e,Y,p,T,c)
  use mpi
  USE globals, ONLY : flen,jobdir
  USE nozfull_data, ONLY: BL_flag 
 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: rho,u,v,w,e

  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: p,T,c

  INTEGER, PARAMETER :: statsUnit=20
  CHARACTER(LEN=flen) :: statsFile,command
  INTEGER :: ios

  ! construct name of stats file
  !WRITE(statsFile,'(2A)') TRIM(jobdir),'/statistics'


  
  ! See if the BL file exists.. if it doesn't, do nothing.  Otherwise, perturb the BL
  WRITE(statsfile,'(2A)') TRIM(jobdir),'/BL.dat'
  OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='OLD',IOSTAT=ios)
  CLOSE(statsUnit)

  ! ios=29 for non-existant file
  ! ios=0 for existant file, turn on the BL flag and delete the file
  IF(ios==0) THEN
      BL_flag = .TRUE.
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      WRITE(command,'(2A)') 'rm -f ', TRIM(statsFile)
      IF(xyzcom_id==0) CALL SYSTEM(command)
  END IF

  
END SUBROUTINE prob_stats

SUBROUTINE prob_source(RHS,simtimeL)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT) :: RHS
  DOUBLE PRECISION, INTENT(IN) :: simtimeL
END SUBROUTINE prob_source

! -----------------------------------------------------------------------------------
! prob_plots
! Problem specific function to write x-y plotfile data.  Called just after
! of graphics dump in "viz()" of file "viz.f".
! -----------------------------------------------------------------------------------
SUBROUTINE prob_plots(plotdir)
  USE mpi
  USE globals, ONLY: flen,ax,ay,az,ix,iy,iz,rho,p,u,v,w,Y
  USE inputs, ONLY: nx,ny,nz,z1,dz,gfilter
  USE prob_interface, ONLY : jobname
  USE interfaces, ONLY : SUMproc,subsum3xy,subsum3yz,subsum3xz,filter
  USE nozfull_data
  USE metrics, ONLY: x_c,y_c,z_c
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: plotdir
  CHARACTER(LEN=flen) :: plotFile
  INTEGER, PARAMETER :: plotUnit=21
  INTEGER :: k
  DOUBLE PRECISION :: x_a,x_b,y_a,y_b,z_a,z_b,core,x0
  DOUBLE PRECISION, DIMENSION(ny) :: ubar,vbar,wbar,rhobar,pbar,ybar
  DOUBLE PRECISION, DIMENSION(nx) :: p_dwn,p_up,p_cen,xbar
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: rhop,up,vp,wp,pp
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: tmp0,tmp1,box
  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u2d,v2d,w2d,p2d,rho2d,x2d,y2d

  ! Pressure down length of nozzle
  z_a = 0.0d0 !! dble(nz-1)*dz / 2.0d0 
  z_b = z_a + dble(nz-1)*dz

  where (z_c >= z_a .and. z_c <= z_a)
      box = 1.0d0
  elsewhere
      box = 0.0d0
  end where

  tmp0 = 0.0d0
  if (iy(1)==1) tmp0(:,1,:) = 1.0d0
  tmp1 = box*tmp0

  p_dwn = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)
  xbar = SUBSUM3YZ(x_c*tmp1)/SUBSUM3YZ(tmp1)

  tmp0 = 0.0d0
  if (iy(ay)==ny) tmp0(:,ay,:) = 1.0d0
  tmp1 = box*tmp0
  p_up = SUBSUM3YZ(p*tmp1)/SUBSUM3YZ(tmp1)

  core = del_BL*2.0D0 
  where (y_c<=-core .and. y_c>=core) box=0.0d0
  p_cen = SUBSUM3YZ(p*box)/SUBSUM3YZ(box)


  ! master cpu writes plot file--------------------------------------------------------------------
  SELECT CASE(xyzcom_id)
  CASE(master)
    WRITE(plotFile,'(2A)') TRIM(plotdir),'/pressure.dat'
    OPEN(UNIT=plotUnit,FILE=TRIM(plotFile),FORM='FORMATTED',STATUS='REPLACE')
    WRITE(plotUnit,*) "%# <1-4> x-loc,p_lower,p_upper,p_cen"
    DO k=1,nx
      WRITE(plotUnit,'(4ES12.4)') xbar(k),p_dwn(k),p_up(k),p_cen(k)
    END DO
    CLOSE(plotUnit)
  END SELECT


END SUBROUTINE prob_plots

! -----------------------------------------------------------------------------------
! prob_bc
! Set the problem-specific boundary values for the state data
! Note that only problem specific BCs need be set.  If the type of the boundary
! is general, such as WALL, SYMM, SLIP, PFLO, PINF, it has already been set.
! This routine is called from "boundary()" in boundary.f
! -----------------------------------------------------------------------------------
SUBROUTINE prob_bc(rho,u,v,w,e,Y)
  USE mpi
  USE globals, ONLY : p,x1proc,y1proc,z1proc,xnproc,ynproc,znproc,molwts,ax,ay,az,iy,ix,iz
  USE globals, ONLY : simtime,step
  USE inputs, ONLY : nx,ny,nz,gamma,gfilter,lfilter
  USE constants
  USE nozfull_data
  USE metrics
  USE interfaces, ONLY : filter !, ran1, SUMproc
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: rho,u,v,w
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: e
  DOUBLE PRECISION, DIMENSION(:,:,:,:), INTENT(INOUT), OPTIONAL :: Y
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: dumF,dumT
  DOUBLE PRECISION, DIMENSION(size(rho,1),size(rho,3))  :: tmpw,blowx,blowz,blow
  DOUBLE PRECISION, DIMENSION(SIZE(u,1),SIZE(u,3)) :: dxdA,dydA,mag,u1,u2,v1,v2	!  A-column of Jacobian Tensor
  !DOUBLE PRECISION, DIMENSION(1) :: randT
  DOUBLE PRECISION :: del
  DOUBLE PRECISION :: filpt,thick
  INTEGER :: i,k


      ! Digital Filtering inflow BC
      CALL DFinflow(rho,u,v,w,e)

      
      IF (xnproc) THEN
         !  Set the outflow BC to be constant
         u(ax,:,:)    = zero
         v(ax,:,:)    = zero
         w(ax,:,:)    = zero
         rho(ax,:,:)  = rho_amb
         e(ax,:,:)    = e_amb
      END IF
      


      !!!  TOP WALL  !!!
      IF (ynproc) THEN

         u(:,ay,:) = zero
         v(:,ay,:) = zero
         w(:,ay,:) = zero

         ! Adiabatic Wall
         !rho(:,ay,:) = (54.0d0*rho(:,ay-1,:)-27.0d0*rho(:,ay-2,:)+6.0d0*rho(:,ay-3,:) )/33.0d0
         !e(:,ay,:) = (54.0d0*e(:,ay-1,:)-27.0d0*e(:,ay-2,:)+6.0d0*e(:,ay-3,:) )/33.0d0 
         rho(:,ay,:) = rho(:,ay-1,:) 
         e(:,ay,:)   = e(:,ay-1,:) 

      END IF

      !!!   BOTTOM WALL   !!!
      IF (y1proc) THEN

         u(:,1,:) = zero
         v(:,1,:) = zero
         w(:,1,:) = zero
         
         ! Adiabatic Wall 
         !rho(:,1,:) = (54.0d0*rho(:,2,:)-27.0d0*rho(:,3,:)+6.0d0*rho(:,4,:) )/33.0d0
         !e(:,1,:) = (54.0d0*e(:,2,:)-27.0d0*e(:,3,:)+6.0d0*e(:,4,:) )/33.0d0 
         rho(:,1,:) = rho(:,2,:) 
         e(:,1,:)   = e(:,2,:) 

      END IF


      !!!   SIDE WALL !!!
      !IF (z1proc) THEN

!         u(:,:,1) = zero
!         v(:,:,1) = zero
!         w(:,:,1) = zero
!         rho(:,:,1) = rho_in !rho(:,:,az-1) 
!         e(:,:,1)   = e_in !e(:,:,az-1) 

!      END IF


      !!!   SIDE WALL !!!
!      IF (znproc) THEN

!         u(:,:,az) = zero
!         v(:,:,az) = zero
!         w(:,:,az) = zero
!         rho(:,:,az) = rho_in !rho(:,:,az-1) 
!         e(:,:,az)   = e_in !e(:,:,az-1) 

!      END IF

      ! Gradually apply the exit boundary conditions  
      filpt = dble(nx-2)
      thick = 3.0d0
      DO i=1,ax
         dumT(i,:,:)=(one+tanh((dble(ix(i))-filpt)/thick))/two
      END DO
      
      !!!  OUT-FLOW  !!!
      IF (xnproc) THEN 
         ! Force back pressure to remain ambient and let density float from NSCBC
         rho = rho + dumT * ( rho_amb - rho )                         
         e = e + dumT * ( e_amb - e )                         
         
      END IF  
      
      ! Gussian Filter for last N points in x-direction (A)
      CALL filter(gfilter,u,dumF)
      u = u + dumT*(dumF-u) 
      CALL filter(gfilter,v,dumF)
      v = v + dumT*(dumF-v)
      CALL filter(gfilter,w,dumF)
      w = w + dumT*(dumF-w)
      CALL filter(gfilter,e,dumF)
      e = e + dumT*(dumF-e)
      CALL filter(gfilter,rho,dumF)
      rho = rho + dumT*(dumF-rho)
      
      IF(nz==1) w = zero

      ! Filter the inflow to smooth stuff out
      !filpt = dble(5)
      !thick = 3.0d0
      !DO i=1,ax
      !   dumT(i,:,:)=(one+tanh((dble(ix(i))-filpt)/thick))/two
      !END DO

      !dumT = 1.D0
      !u = zero
      !v = zero
      !w = zero
      !e = e_amb
      !rho = rho_amb
      ! Gussian Filter for last N points in x-direction (A)
      !CALL filter(gfilter,u,dumF)
      !u = u + dumT*(dumF-u) 
      !CALL filter(gfilter,v,dumF)
      !v = v + dumT*(dumF-v)
      !CALL filter(gfilter,w,dumF)
      !w = w + dumT*(dumF-w)
      !CALL filter(gfilter,e,dumF)
      !e = e + dumT*(dumF-e)
      !CALL filter(gfilter,rho,dumF)
      !rho = rho + dumT*(dumF-rho)



  
END SUBROUTINE prob_bc



! -----------------------------------------------------------------------------------
! prob_acceleration
! Set the body force source term for a compressible calculation
! called by "acceleration" in "acceleration.f"
! -----------------------------------------------------------------------------------
SUBROUTINE prob_acceleration(simtime,gx,gy,gz)
  USE mpi
  USE prob_interface, ONLY : jobname
  USE inputs, ONLY : compressible,accelx,accely,accelz
  USE constants, ONLY : zero
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: gx,gy,gz

  gx = accelx
  gy = accely
  gz = accelz


END SUBROUTINE prob_acceleration


SUBROUTINE prob_esource()

! Do something

END SUBROUTINE prob_esource


! ------------------------------------------------------------------------------
! Set the body force, due to frame acceleration, for incompressible cases.
! ------------------------------------------------------------------------------
 SUBROUTINE prob_Iacceleration(simtime,g)
  USE inputs, ONLY : accelz
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: simtime
  DOUBLE PRECISION, INTENT(OUT) :: g
  g = accelz
 END SUBROUTINE prob_Iacceleration
 
 
 ! ------------------------------------------------------------------------------
! Assign species enthalpies compatible with prob_eos. These are necessary in
! order to compute the enthalpy diffusion term in the energy equation.
! ------------------------------------------------------------------------------
 FUNCTION prob_enthalpies()
  USE constants
  USE inputs, ONLY : ns
  USE globals, ONLY : ax,ay,az,p
  USE leosvars, ONLY : part_r,part_e
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(ax,ay,az,ns) :: prob_enthalpies
  INTEGER :: i
     DO i=1,ns
       prob_enthalpies(:,:,:,i) = part_e(:,:,:,i) + p/part_r(:,:,:,i) ! assumes all species are at the same pressure
     END DO
 END FUNCTION prob_enthalpies

! ------------------------------------------------------------------------------
! Assign problem-specific viscosity, thermal conductivity,
! species diffusivities and magnetic diffusivity.
! ------------------------------------------------------------------------------ 
 SUBROUTINE prob_properties()
  USE constants, ONLY : zero,one,half,three,Runiv
  USE inputs, ONLY : ns,viscous,conductive,diffusive,resistive,diffusive,gamma
  USE globals, ONLY : materials,T
  USE globals, ONLY : molwts
    USE globals, ONLY : mu,bulk,ktc,Diff,eta ! *** OUTPUTS ***
  USE leosvars, ONLY : part_r
  USE nozfull_data
  IMPLICIT NONE
  DOUBLE PRECISION :: T_0,ST,R,cp,cv

  ! Specific Gas Constant
  R = Runiv/molwts(1)
      
  ! Sutherland's Law for Viscosity
  T_0 = 273.15D0  ! Kelvin reference temperature
  ST = 110.4D0    ! Sutherland temperature
  mu = mu_0 * (T/T_0)**(three*half) * (T_0 + ST) / (T + ST)

  ! Bulk viscosity
  bulk = zero
  
  ! Thermal diffusion
  cv = R/(gamma-1.0D0)
  cp = gamma*cv
  ktc = cp * mu / Pr

  IF (ASSOCIATED(Diff)) THEN
    Diff = zero
  END IF
  IF (ASSOCIATED(eta)) THEN
    eta = zero
  END IF
 END SUBROUTINE prob_properties

!-------------------------------------------------------------------------------
! Assign problem-specific opacities. This routine is only called if
! eostype='TABLE' and opacities=2.
! ------------------------------------------------------------------------------ 
 SUBROUTINE prob_opacities()
  USE constants, ONLY : zero,one
  USE inputs, ONLY : ns,radiative
  USE globals, ONLY : materials,rho,p,T,Y
  USE globals, ONLY : atomic_weight,atomic_number,atomic_number23,atomic_number2
  USE globals, ONLY : kP,kR ! *** OUTPUTS ***
  USE leosvars, ONLY : part_r
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Y,1),SIZE(Y,2),SIZE(Y,3)) :: tmp,dum
  DOUBLE PRECISION, DIMENSION(SIZE(Y,1),SIZE(Y,2),SIZE(Y,3),SIZE(Y,4)) :: zeff
  INTEGER :: i,j,k,n
  IF (ASSOCIATED(kP)) THEN
    kP = zero
  END IF
  IF (ASSOCIATED(kR)) THEN
    kR = zero
  END IF
 END SUBROUTINE prob_opacities


!-------------------------------------------------------------------------------
! Assign problem-specific geometry. This routine is only called if
! upon initialization only
! ------------------------------------------------------------------------------ 
SUBROUTINE prob_geom()
 USE constants
 USE metrics, ONLY: x_c,y_c,z_c
 USE nozfull_data
 USE mpi
 USE constants
 USE inputs
 USE globals
 USE interfaces, ONLY : filter
 IMPLICIT NONE
 DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: tmp,dum
 INTEGER :: i,j,k,nz_tmp
 INTEGER :: re_unit,ii

  re_unit = 19
  ALLOCATE(tmp(256,ny,nz))
  ALLOCATE(dum(ax,ay,az))

  !  Put file read routine here 
  ! read X,Y and Z points
  ! Get X
  WRITE(gridfile,'(2A)') TRIM(gridpath),'/X.grd'
  OPEN(UNIT=re_unit,FILE=TRIM(gridfile),FORM='UNFORMATTED',STATUS='OLD')
  READ(re_unit) tmp
  x_c(:,:,:) = tmp(ix,iy,iz)
  CLOSE(re_unit)

  ! Get Y
  WRITE(gridfile,'(2A)') TRIM(gridpath),'/Y.grd'
  OPEN(UNIT=re_unit,FILE=TRIM(gridfile),FORM='UNFORMATTED',STATUS='OLD')
  READ(re_unit) tmp
  y_c(:,:,:) = tmp(ix,iy,iz)
  CLOSE(re_unit)

  ! Get Z
  WRITE(gridfile,'(2A)') TRIM(gridpath),'/Z.grd'
  OPEN(UNIT=re_unit,FILE=TRIM(gridfile),FORM='UNFORMATTED',STATUS='OLD')
  READ(re_unit) tmp
  z_c(:,:,:) = tmp(ix,iy,iz)
  CLOSE(re_unit)

  DEALLOCATE(tmp)

  !  For Nozzle-grid convert mm -> cm
  !  Z is left out because it is specified through dz
  x_c = x_c*1.0d-1
  y_c = y_c*1.0d-1
  z_c = z_c*1.0d-1

  !! Make sure everyone waits till script if finished 
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  !! Filter the grid...
  DO ii=1,3
  dum = x_c
  CALL filter(lfilter,dum,x_c)
  dum = y_c
  CALL filter(lfilter,dum,y_c)
  dum = z_c
  CALL filter(lfilter,dum,z_c)
  END DO

  !! Get the metric terms here.
  CALL get_jacobian()

END SUBROUTINE



!-------------------------------------------------------------------------------
! Assign problem-specific geometry. This routine is only called if
! upon initialization only
! ------------------------------------------------------------------------------ 
SUBROUTINE prob_geom_2d()
 USE constants
 USE metrics, ONLY: x_c,y_c,z_c
 USE nozfull_data
 USE mpi
 USE constants
 USE inputs
 USE globals
 USE interfaces, ONLY : filter
 IMPLICIT NONE
 DOUBLE PRECISION, DIMENSION(ax,ay,az) :: tophat,tmp
 DOUBLE PRECISION, DIMENSION(nx,ny) :: tmpx,tmpy
 INTEGER :: i,j,k,nz_tmp
 INTEGER :: re_unit


  re_unit = 19
  
  
  !  Put file read routine here 
  re_unit = 35
  OPEN(UNIT=re_unit,FILE=TRIM(gridfile),STATUS='OLD')
  DO J=1,ny
      DO I = 1,nx
         READ(re_unit,*) tmpx(I,J),tmpy(I,J)
      END DO
  END DO
   
  x_c(:,:,1) = tmpx(ix,iy)
  y_c(:,:,1) = tmpy(ix,iy)
  CLOSE(re_unit)

  !  Z-direction extrusion
  DO k=1,az
      z_c(:,:,k) = zcom_id*az*dz + dz*dble(k-1)
      x_c(:,:,k) = x_c(:,:,1)
      y_c(:,:,k) = y_c(:,:,1)
  END DO

  !  For Nozzle-grid convert mm -> cm
  !  Z is left out because it is specified through dz
  x_c = x_c*1.0d-1
  y_c = y_c*1.0d-1


  !! Make sure everyone waits till script if finished 
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  ! Make sure Jacobian is 2d... extruded mesh
  nz_tmp = nz
  nz = 1
  CALL get_jacobian()
  nz = nz_tmp

END SUBROUTINE





FUNCTION sech(x)
 DOUBLE PRECISION :: x,sech
      sech = 1.0d0 / cosh(x)
END FUNCTION 


SUBROUTINE interpZv(X1,X2,Y1,Y2,n1,n2,nZ,var)
  !! Assumes there is no variation in x1,x2,y1,y2 int the z (2nd) index direction
  IMPLICIT NONE
  INTEGER :: N1,N2,nz,var
  DOUBLE PRECISION, DIMENSION(N1),    INTENT(IN)  :: X1
  DOUBLE PRECISION, DIMENSION(var,N1,nz), INTENT(IN)  :: Y1
  DOUBLE PRECISION, DIMENSION(N2),    INTENT(IN)  :: X2
  DOUBLE PRECISION, DIMENSION(var,N2,nz), INTENT(OUT) :: Y2
  INTEGER :: j,jm,i
  
      j = 1
      jm= 0 
      DO i=1,n2
         DO WHILE(1==1)
            IF( x2(i) == x1(j) ) THEN ! Are we dead on ?
               y2(:,i,:) = y1(:,j,:)
               EXIT
            ELSEIF( j == n1) THEN ! Were at the end and still havent found a bounds.... just give it last point
               y2(:,i,:) = y1(:,n1,:)
               EXIT
            ELSEIF( x1(j) < x2(i) ) THEN ! Are we too small
               j = j + 1
               jm= j - 1
            ELSEIF( x1(j) > x2(i) .and. x1(jm) < x2(i)  ) THEN ! In range... interpolate
               jm = j - 1
               y2(:,i,:) = y1(:,jm,:) + ( y1(:,j,:) - y1(:,jm,:))/(x1(j)-x1(jm))*(x2(i)-x1(jm))
               EXIT
            ELSE                ! Too far right... move left one
               j = j - 1
               jm= j - 1
            END IF
         END DO
      END DO


END SUBROUTINE


SUBROUTINE interp1D(X1,X2,Y1,Y2,n1,n2)
  !! Assumes there is no variation in x1,x2,y1,y2 int the z (2nd) index direction
  IMPLICIT NONE
  INTEGER :: N1,N2,nz,var
  DOUBLE PRECISION, DIMENSION(N1),    INTENT(IN)  :: X1
  DOUBLE PRECISION, DIMENSION(N1), INTENT(IN)  :: Y1
  DOUBLE PRECISION, DIMENSION(N2),    INTENT(IN)  :: X2
  DOUBLE PRECISION, DIMENSION(N2), INTENT(OUT) :: Y2
  INTEGER :: j,jm,i
  
      j = 1
      jm= 0 
      DO i=1,n2
         DO WHILE(1==1)
            IF( x2(i) == x1(j) ) THEN ! Are we dead on ?
               y2(i) = y1(j)
               EXIT
            ELSEIF( j == n1) THEN ! Were at the end and still havent found a bounds.... just give it last point
               y2(i) = y1(n1)
               EXIT
            ELSEIF( x1(j) < x2(i) ) THEN ! Are we too small
               j = j + 1
               jm= j - 1
            ELSEIF( x1(j) > x2(i) .and. x1(jm) < x2(i)  ) THEN ! In range... interpolate
               jm = j - 1
               y2(i) = y1(jm) + ( y1(j) - y1(jm))/(x1(j)-x1(jm))*(x2(i)-x1(jm))
               EXIT
            ELSE                ! Too far right... move left one
               j = j - 1
               jm= j - 1
            END IF
         END DO
      END DO


END SUBROUTINE



SUBROUTINE random_BL(u,v,w)
      USE mpi
      USE inputs, ONLY: nx,ny,nz,gfilter,lfilter
      USE globals, ONLY: ax,ay,az,ix,iy,iz,y1proc,ynproc,x1proc,xnproc
      USE metrics, ONLY : x_c,y_c
      USE constants, ONLY: zero,half,one,two
      USE interfaces, ONLY: restart,filter,filtery,set_time_seed,ran1,SUBSUM3YZ
      USE nozfull_data
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(ax,ay,az),INTENT(INOUT) :: u,v,w
      DOUBLE PRECISION, DIMENSION(ax,ay,az) :: wU,wL,wN,ywall,weight
      DOUBLE PRECISION, DIMENSION(ax,1,1) :: tmp
      DOUBLE PRECISION, DIMENSION(ay,az) :: up,vp,wp
      DOUBLE PRECISION, DIMENSION(ny) :: y_U,y_L
      DOUBLE PRECISION, DIMENSION(size(u,1),size(u,2),size(u,3))  :: tmpy,tmpx,tmpu
      INTEGER :: in_unit,i,j,k,re_unit,pt
      DOUBLE PRECISION :: dumy,bL,dum1,dum2,shift,umag,BB,f1w,f2w,Uinf,dt_step,thick
      DOUBLE PRECISION :: plane,xBL
      INTEGER :: nxp,nyp,wall
      CHARACTER(LEN=30) :: comments,filename

      !! Use this to get smooth profile at wall
      IF(y1proc) u(:,1,:)  = zero
      IF(ynproc) u(:,ay,:) = zero
      IF(xyzcom_id == 0) print*,'Perturbing BL region'
      
      !DO i=1,10
      !   CALL filtery(gfilter,u,tmpx)
      !   CALL filtery(gfilter,tmpx,u)
      !END DO

      CALL set_time_seed(1)
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      

      !! BJO -- Make the j index dependent on the BL height
      !!  also make BL length only up to the nozzle or ....
      !! Perturb the present flow... dont filter the flow... this smoothes things and will slow transition.
      !
      DO i=1,ax
         DO j=1,ay
            DO k=1,az
               xBL = half*( one - tanh( dble(ix(i)-rcy_pt_g)/ 4.0d0 ) )
               wL(i,j,k) = xBL * half*( one - tanh( dble(iy(j)-30) / 4.0d0 ))
               wU(i,j,k) = xBL * half*( one + tanh( dble(iy(j)-( ny-30) ) / 4.0d0 ))
            END DO 
         END DO
      END DO

      weight = wL + wU

      weight = weight !* .1d0

      f1w = .3d-1    ! High frquency
      f2w = .4d0    ! Low  frequency
      
      tmpu = u
      
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      u = u * (one + weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      u = u * (one + weight * f2w*tmpx) ! /sqrt(tmpy))
      
      
      
      f1w = .1d-1                ! High frquency
      f2w = .2d0                ! Low  frequency

      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      v = v +  tmpu* (weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      v = v + tmpu* (weight * f2w*tmpx) ! /sqrt(tmpy))


      f1w = .1d-1                ! High frquency
      f2w = .2d0                ! Low  frequency

      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(lfilter,tmpy,tmpx)   ! High Wave Number filter
      w = w + tmpu* ( weight * f1w*tmpx) ! /sqrt(tmpy))
      CALL ran1(tmpy)
      tmpy = two*(tmpy-half)           ! Make RN from (-1,1)
      CALL filter(gfilter,tmpy,tmpx)   ! Low Pass filter
      w = w +  tmpu* ( weight * f2w*tmpx) ! /sqrt(tmpy))

      IF(y1proc) THEN
         u(:,1,:) = zero
         v(:,1,:) = zero
         w(:,1,:) = zero
      END IF

      IF(ynproc) THEN
         u(:,ay,:) = zero
         v(:,ay,:) = zero
         w(:,ay,:) = zero
      END IF

      
END SUBROUTINE random_BL







SUBROUTINE get_2d_z(var,plane)
      USE globals, ONLY: ax,ay,az,ix,iy,iz
      USE inputs, ONLY: nx,ny,nz
      USE interfaces, ONLY: SUBSUM3YZ
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(IN) :: var
      DOUBLE PRECISION, DIMENSION(nx,ny) ,INTENT(OUT)   :: plane
      
      DOUBLE PRECISION, DIMENSION(ax,1,az) :: xlineL
      DOUBLE PRECISION, DIMENSION(nx) :: xlineG
      INTEGER :: i,ii

      DO i=1,ny
         ii = MOD( (i-1), ay ) + 1
         
         xlineL = 0.0D0
         IF( iy(ii) == i ) xlineL(:,1,:) = var(:,ii,:)
           
         xlineG = SUBSUM3YZ(xlineL)
         plane(:,i) = xlineG
      END DO

      plane = plane / dble(nz)

END SUBROUTINE



!! Routines for turbulent inflow via digital filtering technique.
SUBROUTINE setup_DFinflow
  USE mpi, ONLY: xyzcom_id,master
  USE inputs, ONLY: ny
  USE globals, ONLY: ay,az,ix,iy,iz,simtime,x1proc,jobdir,dt
  USE constants, ONLY: zero,one,two
  USE interfaces, ONLY: SUBSUM3XZ,newdir
  USE metrics, ONLY: y_c
  USE nozfull_data !, ONLY: T_in,P_in,rho_in,U_in,Mach,A11,A22,A33,A12,Umean,Tmean
  !USE nozfull_data, ONLY: RHO_u,RHO_v,RHO_w,UIx,UIy,UIz,VIx,VIy,VIz,WIx,WIy,WIz
  !USE nozfull_data, ONLY: Nbuff,del_star,del_BL,y_r,y_G,buI,bvI,bwI,buO,bvO,bwO
  !USE nozfull_data, ONLY: simtime_old,tauX

  IMPLICIT NONE
  DOUBLE PRECISION :: line,in_plane
  DOUBLE PRECISION, DIMENSION(ay,az) :: UUmean,VVmean,WWmean,UVmean
  DOUBLE PRECISION, DIMENSION(1,ay,1) :: tmp
  DOUBLE PRECISION, DIMENSION(ny) :: Uprof
  DOUBLE PRECISION, DIMENSION(ny/2) :: Utmp
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yData,vData     
  INTEGER :: runit=43
  INTEGER :: ndata,i,isync
  CHARACTER(len=80) :: uaveFile

      ! Allocate the mean profiles
      ALLOCATE(Umean(ay,az))
      ALLOCATE(Tmean(ay,az))
      ALLOCATE(A11(ay,az))
      ALLOCATE(A22(ay,az))
      ALLOCATE(A33(ay,az))
      ALLOCATE(A12(ay,az))

      ALLOCATE(RHO_u(ay,az))
      ALLOCATE(RHO_v(ay,az))
      ALLOCATE(RHO_w(ay,az))
      
      ALLOCATE(y_G(ny))
      ALLOCATE(y_r(ay))



      in_plane = zero
      line = zero
      IF (ix(1) == 1) in_plane = one
      IF (iz(1) == 1) line = one
      tmp(1,:,1) = y_c(1,:,1)*in_plane*line
      y_G = SUBSUM3XZ(tmp)
      y_r = y_G(1:ay)           ! Only take the profile on this proc.
      
      IF ( iy(1) < ny/2) THEN
         y_r = y_r - y_G(1)     ! Bottom boundary layer
      ELSE
         y_r = y_G(ny) - y_r    ! Top boundary layer
      END IF

      y_G = y_G - y_G(1)  ! 1-ny ascending.  Do one side then mirror mean
      y_G = y_G / del_BL   ! Non-dimensionalize by BL thickness... for interp only
      ndata = 384/2       ! Length of the profile data array to read in
      ALLOCATE(yData(ndata))
      ALLOCATE(vData(ndata))

      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/Uprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')

      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         Umean(:,i) = U_in * Uprof(iy)
      END DO

      ! Calculate the Displacement thickness
      del_star = zero
      DO i=1,ny/2
         del_star = del_star + (one - Uprof(i))* (y_G(i+1)-y_G(i))
      END DO
      del_star = del_star * del_BL  ! The U prof was non-dim by del_BL


      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/Tprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')
      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         Tmean(:,i) = T_in * Uprof(iy)
      END DO

      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/UUprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')
      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         UUmean(:,i) = U_in**two * Uprof(iy)
      END DO

      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/VVprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')
      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         VVmean(:,i) = U_in**two * Uprof(iy)
      END DO

      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/WWprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')
      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         WWmean(:,i) = U_in**two * Uprof(iy)
      END DO

      ! Set up some profiles here
      OPEN(UNIT=runit,FILE='FINE/UVprof.dat',STATUS='OLD', & 
      & FORM='FORMATTED')
      DO i=1,ndata
         READ(runit,*) vData(i),yData(i)
      END DO
      CLOSE(runit)
      CALL interp1D(yData,y_G(1:ny/2),vData,Utmp,ndata,ny/2)
      Uprof(1:ny/2) = Utmp
      Uprof(ny:ny/2+1:-1) = Utmp
      DO i=1,az
         UVmean(:,i) = U_in**two * Uprof(iy)
      END DO

      ! Get the Reynolds Stress tensor
      A11 = sqrt(UUmean)
      A12 = UVmean / A11
      WHERE (A11 == zero) A12 = zero  ! Treat zero 
      A22 = sqrt( VVmean - A12**two )
      A33 = sqrt(WWmean)

      ! Non-dimensionalize the local wall normal length scale for filter size selection
      y_r = y_r / del_star


      ! Setup the 2d filters ( 6 total.  Inner (u,v,w) and outer (u,v,w) )
      UIx = 10
      UIy = (/ 20, 35 /)
      UIz = 20

      VIx = 4
      VIy = (/ 25, 45 /)
      VIz = 20
      
      WIx = 4
      WIy = (/ 15, 20 /)
      WIz = 30

      Nbuff = 45 * 2

      simtime_old = simtime
      tauX = DBLE(UIx)*del_star / U_in

      ! Allocate the filter coefficients
      ALLOCATE(buI(2*UIy(1)+1,2*UIz+1))
      ALLOCATE(bvI(2*VIy(1)+1,2*VIz+1))
      ALLOCATE(bwI(2*WIy(1)+1,2*WIz+1))

      ALLOCATE(buO(2*UIy(2)+1,2*UIz+1))
      ALLOCATE(bvO(2*VIy(2)+1,2*VIz+1))
      ALLOCATE(bwO(2*WIy(2)+1,2*WIz+1))


      ! Inner Filters
      CALL setFiltCoeffs(UIy(1),UIz,buI)
      CALL setFiltCoeffs(VIy(1),VIz,bvI)
      CALL setFiltCoeffs(WIy(1),WIz,bwI)

      ! Outer Filters
      CALL setFiltCoeffs(UIy(2),UIz,buO)
      CALL setFiltCoeffs(VIy(2),VIz,bvO)
      CALL setFiltCoeffs(WIy(2),WIz,bwO)

      ! Set up directory 
      WRITE (uaveDir,'(2A)') TRIM(jobdir),'/Dfil'
      CALL newdir(LEN_TRIM(uaveDir),TRIM(uaveDir),isync)
         

      ! Initialize the temporal averages here.  They will be 
      ! over written on restart read.
      RHO_u = zero
      RHO_v = zero
      RHO_w = zero

      ! Read File for time averages when NOT at t0=0
      IF (simtime .GT. dt .and. init_stats .eq. .FALSE.) THEN
         IF(xyzcom_id == master) PRINT*,'Reading old DF Files'
            
         ! READ FILE only if on inlet BC
         IF (x1proc) THEN
            WRITE (uaveFile,'(2A,I6.6)') TRIM(uaveDir),'/p',xyzcom_id
            runit = 37
            OPEN(UNIT=runit,FILE=TRIM(uaveFile),FORM='UNFORMATTED', &
            & STATUS='UNKNOWN')
            READ(runit) RHO_u
            READ(runit) RHO_v
            READ(runit) RHO_w
            CLOSE(runit)
         END IF
      END IF

      
END SUBROUTINE setup_DFinflow
      
SUBROUTINE get_rands(ny,nz,Nbuff,rands)
  USE mpi, ONLY: xyzcom_id,master
  USE interfaces, ONLY: SUMproc,ran1
  USE ran_state
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ny,nz,Nbuff
  DOUBLE PRECISION, DIMENSION(4,ny+Nbuff,nz+Nbuff),INTENT(OUT) :: rands

  INTEGER :: time_seed 
  INTEGER, DIMENSION(8) :: dtval
  DOUBLE PRECISION :: dubSeed

  ! Master CPU gets the seed for this sequence
  time_seed = 0
  IF (xyzcom_id == master) THEN
      CALL system_clock(count=time_seed)    
      time_seed = 2*(xyzcom_id+1)*abs(time_seed/2)+1
  END IF
  
  ! Send this seed to all procs
  dubSeed = time_seed  
  dubSeed = SUMproc(dubSeed)
  time_seed = INT(dubSeed)

  ! Get the RNs... same for all procs
  CALL ran_seed(sequence=time_seed)  
  CALL ran1(rands)

END SUBROUTINE get_rands


SUBROUTINE DFinflow(rho,u,v,w,e)
  USE mpi
  USE globals, ONLY: x1proc,simtime,ax,ay,az,dump
  USE inputs, ONLY: gamma,nx,ny,nz,gfilter
  USE constants, ONLY: zero,one,two,pi,half
  USE nozfull_data 
  USE interfaces, ONLY: filter,gaufily,gaufilz,MAXVAL2D
  

  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(INOUT) :: rho,u,v,w,e
  
  DOUBLE PRECISION, DIMENSION(ax,ay,az) :: T,P
  DOUBLE PRECISION, DIMENSION(ay,az) :: uinlet,vinlet,winlet,rhoprime,Tprime,MaSq
  DOUBLE PRECISION, DIMENSION(ay,az) :: vU,vV,vW
  DOUBLE PRECISION, DIMENSION(1,ay,az) :: filtmp,fildum
  DOUBLE PRECISION, DIMENSION(4,ny+Nbuff,nz+Nbuff) :: rands,rtmp
  DOUBLE PRECISION :: gm1,EXPt,t_nm1,t_n,dt_n
  INTEGER :: nx_tmp,runit
  CHARACTER(len=80) :: uaveFile

  ! Thermo variables
  gm1 = gamma - one


  ! Get random numbers, need 2 sets, of size (ny+buff,nz+buff)... buff is for the filter width
  CALL get_rands(ny,nz,Nbuff,rands)


  ! Make them have normal distribution (Box-Mueller theorem)
  rtmp(1:2,:,:) = sqrt( -two*LOG(rands((/1,3/),:,:))) * cos(two*pi*rands((/2,4/),:,:))
  rtmp(3:4,:,:) = sqrt( -two*LOG(rands((/1,3/),:,:))) * sin(two*pi*rands((/2,4/),:,:))
  rands = rtmp

  ! Filter them here
  CALL filtRands(UIz,UIy(1),UIy(2),Nbuff,buI,buO,rands(1,:,:),vU)
  CALL filtRands(VIz,VIy(1),VIy(2),Nbuff,bvI,bvO,rands(2,:,:),vV)
  CALL filtRands(WIz,WIy(1),WIy(2),Nbuff,bwI,bwO,rands(3,:,:),vW)

  
  ! Check to see if dump has incremented
  ! If it has, write new Uave file before they change for this next time step
  IF (res_dump .NE. dump) THEN
      ! WRITE FILE
      IF (x1proc) THEN
         WRITE (uaveFile,'(2A,I6.6)') TRIM(uaveDir),'/p',xyzcom_id
         runit = 17
         OPEN(UNIT=runit,FILE=TRIM(uaveFile),FORM='UNFORMATTED', &
         & STATUS='REPLACE')
         WRITE(runit) RHO_u
         WRITE(runit) RHO_v
         WRITE(runit) RHO_w
         CLOSE(runit)
      END IF
  END IF
  res_dump = dump



  ! Get the updated rho_k
  ! Time avergaging coefficients and quantities/fluctuations
  ! Need to write a restart file with averages stored
  t_nm1 = simtime_old 
  t_n = simtime
  simtime_old = simtime
  dt_n = t_n - t_nm1
  EXPt = exp( -pi*dt_n/tauX )
  RHO_u = RHO_u * sqrt( EXPt ) + vU * sqrt( one - EXPt ) 
  RHO_v = RHO_v * sqrt( EXPt ) + vV * sqrt( one - EXPt ) 
  RHO_w = RHO_w * sqrt( EXPt ) + vW * sqrt( one - EXPt ) 

  ! Add perturbations to mean with given 2 point correlations
  uinlet = Umean +  A11 * RHO_u
  vinlet =          A12 * RHO_u + A22 * RHO_v  
  winlet =          A33 * RHO_w
 
  ! Get the temperature perturbation using strong Reynolds Analogy (SRA)
  MaSq = Mach**two !* Umean**two / Tmean
  Tprime = Tmean*( -gm1*MaSq * (uinlet-Umean) / U_in )

  ! Pressure is constant
  IF (x1proc) THEN
      u(1,:,:) = uinlet
      v(1,:,:) = vinlet
      w(1,:,:) = winlet
      T(1,:,:) = Tmean + Tprime
      P(1,:,:) = P_in
      rho(1,:,:) = P_in/(Rgas*T(1,:,:))
      e(1,:,:) = (P_in/(gamma-one))/rho(1,:,:)
  END IF

END SUBROUTINE DFinflow
      



SUBROUTINE filtRands(Nspan,Ni,No,Nbuff,bmnI,bmnO,rands,vfilt)
  USE inputs, ONLY: nx,ny,nz
  USE globals, ONLY: iy,iz,ay,az
  USE nozfull_data, ONLY: y_r
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: Nspan,Ni,No,Nbuff
  DOUBLE PRECISION, INTENT(IN) ::  rands(ny+Nbuff,nz+Nbuff)
  DOUBLE PRECISION, INTENT(IN) ::  bmnI(2*Ni+1,2*Nspan+1), bmnO(2*No+1,2*Nspan+1)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: filt
  DOUBLE PRECISION :: vfilt(ay,az)
  INTEGER :: N1,N2
  INTEGER :: j,k,m,n,mm,nn
  INTEGER :: mF,mG,nF,nG

  N2 = Nspan
  DO j=1,ay

      ! Inner Layer
      IF( y_r(j) < 1.0D0 ) THEN
         N1 = Ni
         IF(.NOT. ALLOCATED(filt)) ALLOCATE(filt(SIZE(bmnI,1),SIZE(bmnI,2) )  )
         filt = bmnI
      ! Outer Layer
      ELSE
         N1 = No
         IF(.NOT. ALLOCATED(filt)) ALLOCATE(filt(SIZE(bmnO,1),SIZE(bmnO,2) ) )
         filt = bmnO
      END IF

      DO k=1,az
         vfilt(j,k) = 0.0D0

         DO m=-N1,N1,1
            mG = (iy(1)-1) + j + (m + N1)
            mF = m + N1 + 1

            DO n=-N2,N2,1
               nG = (iz(1)-1) + k + (n + N2)
               nF = n + N2 + 1
               vfilt(j,k) = vfilt(j,k) + filt(mF,nF)*rands(mG,nG)

            END DO
         END DO
      END DO

  DEALLOCATE(filt)  ! Every new y, compute a new filter
  END DO




END SUBROUTINE filtRands


SUBROUTINE calcDen(N,den)
  USE constants, ONLY: pi,zero,two
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(OUT) :: den
  INTEGER :: i

  den = zero
  DO i=-N,N
      den = den + (exp(-pi*abs(dble(i))/(dble(N)/two)))**two
  END DO

  den = sqrt(den)

END SUBROUTINE calcDen


SUBROUTINE setFiltCoeffs(N1,N2,bij)
  USE constants, ONLY: pi,two
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N1, N2
  DOUBLE PRECISION, DIMENSION(2*N1+1,2*N2+1), INTENT(OUT) :: bij

  INTEGER :: i,j,k1,k2
  DOUBLE PRECISION :: den1,den2,N1d2,N2d2
  
  CALL calcDen(N1,den1)
  CALL calcDen(N2,den2)
  N1d2 = dble(N1)/two
  N2d2 = dble(N2)/two

  DO i=1,2*N1+1
      DO j=1,2*N2+1
         k1 = abs(i-(N1+1))
         k2 = abs(j-(N2+1))
         bij(i,j) = exp(-pi*dble(k1)/N1d2)/den1 * exp(-pi*dble(k2)/N2d2)/den2
      END DO
  END DO

END SUBROUTINE setFiltCoeffs
