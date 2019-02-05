!>Subroutine responsible for guiding the calculation of the fluxes
!@note: Interface 1 is the left hand interface for cell 1
Subroutine InviscidFlux(rLocalSol3D,rLocalFlux3D,rLocalX,rLocalY,rLocalZ,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalLowMachScheme,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalSpecificGasConstant,rLocalCv,rLocalCp,iLocalRecons,iLocalInvOrder)
implicit none
!> Local array sizes
integer, intent(in):: iLocalNcHalo,iLocalNnHalo,iLocalNVar,iLocalNVarComputed,iLocalNnx,iLocalNny,iLocalNnz,iLocalNumberofSpecies
!> Local flag for equation of state
integer, intent(in):: iLocalEquationOfState
!> order of spatial discretisation for inviscid flux
integer, intent(in) :: iLocalInvOrder
!> Local flag for low Mach scheme
integer, intent(in)  :: iLocalLowMachScheme
!> Local solution array
real(8), dimension (iLocalNVar,-iLocalNcHalo+1:iLocalNnx+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNny+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNnz+iLocalNcHalo-1), intent(inout)  :: rLocalSol3D
!> Local 3D flux array
real(8), dimension (iLocalNVar,-iLocalNcHalo+1:iLocalNnx+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNny+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNnz+iLocalNcHalo-1), intent(out)  :: rLocalFlux3D
!> Local X array
real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalX
!> Local Y array
real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalY
!> Local Z array
real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalZ
!> Local gas properties
real(8),dimension (iLocalNumberOfSpecies),intent(in)::rLocalSpecificGasConstant
real(8),dimension (iLocalNumberOfSpecies,2,5),intent(in)::rLocalCv,rLocalCp
!> Flag for reconstruction approach
integer,intent(in)::iLocalRecons
!> Fluxes for each equation
real(8),dimension (4+iLocalNumberOfSpecies)::rLocalFlux
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo)::rStencil
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2)::rLimitedVariables
!> Left and right limited primitive variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2)::rLimitedPrimitiveVariables
!> Left and right limited kinetic energy
real(8),dimension (2)::rLimitedKineticEnergy
! >Metric vector magnitude
real(8) rBeta
!> variable indices
integer i,j,k,v,iLoop
!> Normalised Xi metrics
real(8) rNJDXiDX,rNJDXiDY,rNJDXiDZ
!> Normalised Eta metrics
real(8) rNJDEtDX,rNJDEtDY,rNJDEtDZ
!> Normalised Zeta metrics
real(8) rNJDZeDX,rNJDZeDY,rNJDZeDZ
!> half the stencil size (note that stencils must be an even number of points)
integer iHalfStencilSize
!> One over density
real(8):: rInverseDensity
!> Direction
integer::iDir

if(iLocalRecons.eq.1)then
   !Convert momenta  to primitive variables, i.e. velocity
   do k=-iLocalNcHalo+1,iLocalNnz+iLocalNcHalo-1
      do j=-iLocalNcHalo+1,iLocalNny+iLocalNcHalo-1
         do i=-iLocalNcHalo+1,iLocalNnx+iLocalNcHalo-1
            rInverseDensity=1.d0/rLocalSol3D(iLocalNVar-3,i,j,k)
            do v=1,3
               rLocalSol3D(v,i,j,k)=rLocalSol3D(v,i,j,k)*rInverseDensity
            end do
         end do
      end do
   end do
end if

rLimitedVariables=0.d0
rLimitedPrimitiveVariables=0.d0

!Calculate Xi
idir=1
do k=1,iLocalNnz-1
   do j=1,iLocalNny-1
      do i=1,iLocalNnx
         rLocalFlux=0.d0
         call ConstructXiStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
         ! call Construct_Xi_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rNJDXiDX,rNJDXiDY,rNJDXiDZ,rBeta)
         rNJDXiDX=1.d0
         rNJDXiDY=0.d0
         rNJDXiDZ=0.d0
         call Construct_Cartesian_Xi_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
         call VariableReconstruction(rStencil,rLimitedVariables,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalInvOrder,iHalfStencilSize,iLocalRecons)
         call ComputePCVariables(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant,iLocalRecons)     
         if(iLocalInvOrder.gt.1.and.iLocalLowMachScheme.eq.1)call LowMachCorrection(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalSpecificGasConstant,rLocalCv,rLocalCp)
         call RiemannSolver(rLimitedVariables,rLimitedPrimitiveVariables,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNumberOfSpecies,iLocalEquationOfState,iHalfStencilSize,rLocalSpecificGasConstant,rLocalCv,rLocalCp,rNJDXiDX,rNJDXiDY,rNJDXiDZ,rBeta,rLocalFlux,i,j,k,iDir)
         do v=1,4+iLocalNumberOfSpecies
            !Add the flux to the cell to the left of this interface
            rLocalFlux3D(v,i-1,j,k)=rLocalFlux3D(v,i-1,j,k)+rLocalFlux(v)
            !Subtract the flux from the cell to the right of this interface
            rLocalFlux3D(v,i  ,j,k)=rLocalFlux3D(v,i  ,j,k)-rLocalFlux(v)
         end do
      end do
   end do
end do

!Calculate Eta
idir=2
do k=1,iLocalNnz-1
   do i=1,iLocalNnx-1
      do j=1,iLocalNny
         rLocalFlux=0.d0
         call ConstructEtaStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
         ! call Construct_Et_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rNJDEtDX,rNJDEtDY,rNJDEtDZ,rBeta)
         rNJDEtDX = 0.d0
         rNJDEtDY = 1.d0
         rNJDEtDZ = 0.d0
         call Construct_Cartesian_Et_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
         call VariableReconstruction(rStencil,rLimitedVariables,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalInvOrder,iHalfStencilSize,iLocalRecons)
         call ComputePCVariables(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant,iLocalRecons)
         if(iLocalInvOrder.gt.1.and.iLocalLowMachScheme.eq.1)call LowMachCorrection(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalSpecificGasConstant,rLocalCv,rLocalCp)
         call RiemannSolver(rLimitedVariables,rLimitedPrimitiveVariables,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNumberOfSpecies,iLocalEquationOfState,iHalfStencilSize,rLocalSpecificGasConstant,rLocalCv,rLocalCp,rNJDEtDX,rNJDEtDY,rNJDEtDZ,rBeta,rLocalFlux,i,j,k,iDir)
         do v=1,4+iLocalNumberOfSpecies
            !Add the flux to the cell to the left of this interface
            rLocalFlux3D(v,i,j-1,k)=rLocalFlux3D(v,i,j-1,k)+rLocalFlux(v)
            !Subtract the flux from the cell to the right of this interface
            rLocalFlux3D(v,i,j  ,k)=rLocalFlux3D(v,i,j  ,k)-rLocalFlux(v)
         end do
      end do
   end do
end do

!Calculate Zeta
idir=3
do j=1,iLocalNny-1
   do i=1,iLocalNnx-1
      do k=1,iLocalNnz
         rLocalFlux=0.d0
         call ConstructZetaStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
         ! call Construct_Ze_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rNJDZeDX,rNJDZeDY,rNJDZeDZ,rBeta)
         rNJDZeDX = 0.d0
         rNJDZeDY = 0.d0
         rNJDZeDZ = 1.d0
         call Construct_Cartesian_Ze_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
         call VariableReconstruction(rStencil,rLimitedVariables,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalInvOrder,iHalfStencilSize,iLocalRecons)
         call ComputePCVariables(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant,iLocalRecons)
         if(iLocalInvOrder.gt.1.and.iLocalLowMachScheme.eq.1)call LowMachCorrection(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalSpecificGasConstant,rLocalCv,rLocalCp)
         call RiemannSolver(rLimitedVariables,rLimitedPrimitiveVariables,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNumberOfSpecies,iLocalEquationOfState,iHalfStencilSize,rLocalSpecificGasConstant,rLocalCv,rLocalCp,rNJDZeDX,rNJDZeDY,rNJDZeDZ,rBeta,rLocalFlux,i,j,k,iDir)
         do v=1,4+iLocalNumberOfSpecies
            !Add the flux to the cell to the left of this interface
            rLocalFlux3D(v,i,j,k-1)=rLocalFlux3D(v,i,j,k-1)+rLocalFlux(v)
            !Subtract the flux from the cell to the right of this interface
            rLocalFlux3D(v,i,j,k  )=rLocalFlux3D(v,i,j,k  )-rLocalFlux(v)
         end do
      end do
   end do
end do

if(iLocalRecons.eq.1)then
   !Convert velocities back to momenta
   do k=-iLocalNcHalo+1,iLocalNnz+iLocalNcHalo-1
      do j=-iLocalNcHalo+1,iLocalNny+iLocalNcHalo-1
         do i=-iLocalNcHalo+1,iLocalNnx+iLocalNcHalo-1
            do v=1,3
               rLocalSol3D(v,i,j,k)=rLocalSol3D(v,i,j,k)*rLocalSol3D(iLocalNVar-3,i,j,k)
            end do
         end do
      end do
   end do
end if
end Subroutine InviscidFlux

!>subroutine responsible for constructing the one dimensional stencil in the Xi direction needed for variable reconstruction.
!@note: immediate l and right of the interface is located in rStencil(:,(iLocalNcHalo+1)/2) and rStencil(:,(iLocalNcHalo+1)/2+1)
Subroutine ConstructXiStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
implicit none
!> variable indices
integer i,j,k,v
!> stencil index
integer s
!> half the stencil size (note that stencils must be an even number of points
integer iHalfStencilSize
!> order of spatial discretisation for inviscid flux
integer, intent(in) :: iLocalInvOrder
!>Local array sizes
integer, intent(in):: iLocalNcHalo,iLocalNVar,iLocalNVarComputed,iLocalNnx,iLocalNny,iLocalNnz
!>Local solution array
real(8), dimension (iLocalNVar,-iLocalNcHalo+1:iLocalNnx+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNny+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNnz+iLocalNcHalo-1), intent(in)  :: rLocalSol3D
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo),intent(out)::rStencil
!>position in the stencil
integer iPos

!start of the stencil is at cell (i-iHalfStencilSize) for interface i
select case(iLocalInvOrder)
   case(1)
      iHalfStencilSize=1
   case(2:4)
      iHalfStencilSize=2
   case(5)
      iHalfStencilSize=3
end select

!construct the stencil
do s=1,2*iHalfStencilSize
   iPos=i-1-iHalfStencilSize+s
   do v=1,iLocalNVar!Computed
      !Note here - the (i-1) in rLocalSol3D is to compensate for the fact that s starts from 1
      rStencil(v,s)=rLocalSol3D(v,iPos,j,k)
   end do
end do

end Subroutine ConstructXiStencil

!>subroutine responsible for constructing the one dimensional stencil in the Eta direction needed for variable reconstruction.
!@note: immediate l and right of the interface is located in rStencil(:,(iLocalNcHalo+1)/2) and rStencil(:,(iLocalNcHalo+1)/2+1)
Subroutine ConstructEtaStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
implicit none
!> variable indices
integer i,j,k,v
!> stencil index
integer s
!> half the stencil size (note that stencils must be an even number of points
integer iHalfStencilSize
!> order of spatial discretisation for inviscid flux
integer, intent(in) :: iLocalInvOrder
!>Local array sizes
integer, intent(in):: iLocalNcHalo,iLocalNVar,iLocalNVarComputed,iLocalNnx,iLocalNny,iLocalNnz
!>Local solution array
real(8), dimension (iLocalNVar,-iLocalNcHalo+1:iLocalNnx+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNny+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNnz+iLocalNcHalo-1), intent(in)  :: rLocalSol3D
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo),intent(out)::rStencil
!>position in the stencil
integer jPos

!start of the stencil is at cell (j-iHalfStencilSize) for interface j
select case(iLocalInvOrder)
   case(1)
      iHalfStencilSize=1
   case(2:4)
      iHalfStencilSize=2
   case(5)
      iHalfStencilSize=3
end select

!construct the stencil
do s=1,2*iHalfStencilSize
   jPos=j-1-iHalfStencilSize+s
   do v=1,iLocalNVar!Computed
      !Note here - the (j-1) in rLocalSol3D is to compensate for the fact that s starts from 1
      rStencil(v,s)=rLocalSol3D(v,i,jPos,k)
   end do
end do

end Subroutine ConstructEtaStencil

!>subroutine responsible for constructing the one dimensional stencil in the Zeta direction needed for variable reconstruction.
!@note: immediate l and right of the interface is located in rStencil(:,(iLocalNcHalo+1)/2) and rStencil(:,(iLocalNcHalo+1)/2+1)
Subroutine ConstructZetaStencil(i,j,k,rLocalSol3D,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNnx,iLocalNny,iLocalNnz,iLocalInvOrder,iHalfStencilSize)
implicit none
!> variable indices
integer i,j,k,v
!> stencil index
integer s
!> half the stencil size (note that stencils must be an even number of points
integer iHalfStencilSize
!> order of spatial discretisation for inviscid flux
integer, intent(in) :: iLocalInvOrder
!>Local array sizes
integer, intent(in):: iLocalNcHalo,iLocalNVar,iLocalNVarComputed,iLocalNnx,iLocalNny,iLocalNnz
!>Local solution array
real(8), dimension (iLocalNVar,-iLocalNcHalo+1:iLocalNnx+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNny+iLocalNcHalo-1,-iLocalNcHalo+1:iLocalNnz+iLocalNcHalo-1), intent(in)  :: rLocalSol3D
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo),intent(out)::rStencil
!>position in the stencil
integer kPos

!start of the stencil is at cell (k-iHalfStencilSize) for interface k
select case(iLocalInvOrder)
   case(1)
      iHalfStencilSize=1
   case(2:4)
      iHalfStencilSize=2
   case(5)
      iHalfStencilSize=3
end select

!construct the stencil
do s=1,2*iHalfStencilSize
   kPos=k-1-iHalfStencilSize+s
   do v=1,iLocalNVar!Computed
      !Note here - the (k-1) in rLocalSol3D is to compensate for the fact that s starts from 1
      rStencil(v,s)=rLocalSol3D(v,i,j,kPos)
   end do
end do

end Subroutine ConstructZetaStencil

!>subroutine responsible for reconstructing and limiting (if necessary) the cell interface values (left and right) using the generated stencil
!@note: immediate l and right of the interface is located in rStencil(:,(iLocalNcHalo+1)/2) and rStencil(:,(iLocalNcHalo+1)/2+1)
Subroutine VariableReconstruction(rStencil,rLimitedVariables,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalInvOrder,iHalfStencilSize,iLocalRecons)
implicit none
!>Local array sizes
integer, intent(in):: iLocalNcHalo,iLocalNVar,iLocalNVarComputed
!> order of spatial discretisation for inviscid flux
integer, intent(in) :: iLocalInvOrder
!> Flag for reconstruction approach
integer,intent(in) ::iLocalRecons
!> half the stencil size (note that stencils must be an even number of points
integer,intent(in) :: iHalfStencilSize
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo),intent(in)::rStencil
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(out)::rLimitedVariables
!>variable indices
integer::v
!>Differences
real(8)::rDLL,rDL,rDM,rDR,rDRR,rDLI,rDRI
!>Ratio of Differences
real(8)::rRLL,rRL,rRR,rRRR
!>Limiter parameter
real(8)::rLIMR,rLIML,rBetaR,rBetaL
!>1/30
real(8)::rC130
!>small number
real(8)::rSmall

rSmall=1.d-15

select case(iLocalInvOrder)
   case(1)
      do v=1,iLocalNVarComputed
         !first order reconstruction
         rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
      END DO

      if(iLocalRecons.eq.1)then
         v=iLocalNvar-1
         !first order reconstruction
         rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
      end if

   case(2)
      !Minmod
      do v=1,iLocalNVarComputed
         rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
         rDM=rStencil(v,3)-rStencil(v,2)  !R-L
         rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

         if(rDL.eq.0.d0)then
            rLimL=1.d0 !Maximum value of rLimL=2
         else
            rLIML=MAX(0.d0,MIN(1.d0,rDM/rDL))
         end if

         if(rDR.eq.0.d0)then
            rLimR=1.d0 !Maximum value of rLimR
         else
            rLIMR=MAX(0.d0,MIN(1.d0,rDM /rDR))
         end if

         rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
         rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
      END DO

      if(iLocalRecons.eq.1)then
         v=iLocalNvar-1
         rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
         rDM=rStencil(v,3)-rStencil(v,2)  !R-L
         rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

         if(rDL.eq.0.d0)then
            rLimL=1.d0 !Maximum value of rLimL=2
         else
            rLIML=MAX(0.d0,MIN(1.d0,rDM/rDL))
         end if

         if(rDR.eq.0.d0)then
            rLimR=1.d0 !Maximum value of rLimR
         else
            rLIMR=MAX(0.d0,MIN(1.d0,rDM /rDR))
         end if

         rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
         rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
      end if

   case(3)
      !van Leer
      do v=1,iLocalNVarComputed
         rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
         rDM=rStencil(v,3)-rStencil(v,2)  !R-L
         rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

         if(abs(rDL).lt.rSmall)then
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         else
            rRL =rDM /rDL
            rLIML=2.d0*rRL/(1.d0+rRL)!(rRL +ABS(rRL))/(1.d0+rRL)
            if(rRL.lt.0.d0)rLimL=0.d0
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
         end if

         if(abs(rDR).lt.rSmall)then
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
         else
            rRR =rDM /rDR
            rLimR=2.d0*rRR/(1.d0+rRR)!(rRR +ABS(rRR))/(1.d0+rRR)
            if(rRR.lt.0.d0)rLimR=0.d0
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
         end if

      END DO

      if(iLocalRecons.eq.1)then
         v=iLocalNvar-1
         rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
         rDM=rStencil(v,3)-rStencil(v,2)  !R-L
         rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

         if(abs(rDL).lt.rSmall)then
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         else
            rRL =rDM /rDL
            rLIML=2.d0*rRL/(1.d0+rRL)!(rRL +ABS(rRL))/(1.d0+rRL)
            if(rRL.lt.0.d0)rLimL=0.d0
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
         end if

         if(abs(rDR).lt.rSmall)then
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
         else
            rRR =rDM /rDR
            rLimR=2.d0*rRR/(1.d0+rRR)!(rRR +ABS(rRR))/(1.d0+rRR)
            if(rRR.lt.0.d0)rLimR=0.d0
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
         end if
      end if
   case(4)
      ! Superbee
      do v=1,iLocalNVarComputed
            rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
            rDM=rStencil(v,3)-rStencil(v,2)  !R-L
            rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

            if(rDL.eq.0.d0)then
               rLimL = 1.d0
            else
               rRL = rDM/rDL
               rLimL = max(0.d0,max(min(2.d0,rRL),min(1.d0,2.d0*rRL)))
            end if

            if(rDR.eq.0.d0)then
               rLimR = 1.d0
            else
               rRR = rDM/rDR
               rLimR = max(0.d0,max(min(2.d0,rRR),min(1.d0,2.d0*rRR)))
            end if

            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
         END DO

         if(iLocalRecons.eq.1)then
            v=iLocalNvar-1
            rDL=rStencil(v,2)-rStencil(v,1)  !L-LL
            rDM=rStencil(v,3)-rStencil(v,2)  !R-L
            rDR=rStencil(v,4)-rStencil(v,3)  !RR-R

            if(rDL.eq.0.d0)then
               rLimL = 1.d0
            else
               rRL = rDM/rDL
               rLimL = max(0.d0,max(min(2.d0,rRL),min(1.d0,2.d0*rRL)))
            end if

            if(rDR.eq.0.d0)then
               rLimR = 1.d0
            else
               rRR = rDM/rDR
               rLimR = max(0.d0,max(min(2.d0,rRR),min(1.d0,2.d0*rRR)))
            end if

            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLIML*rDL
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLIMR*rDR
         end if

   case(5)
      !5th order KK

      rC130=1.d0/30.d0

      do v=1,iLocalNVarComputed

         rDLL=rStencil(v,2)-rStencil(v,1) !LL-LLL - i-3/2
         rDL=rStencil(v,3)-rStencil(v,2)  !L-LL - i-1/2
         rDM=rStencil(v,4)-rStencil(v,3)  !R-L - i+1/2
         rDR=rStencil(v,5)-rStencil(v,4)  !RR-R - i+3/2
         rDRR=rStencil(v,6)-rStencil(v,5) !RRR-RR - i+5/2

         if(abs(rDL).lt.rSmall)then !eq.0.d0)then
            !BetaL will be infinite giving rLimL=2 OR 0.d0 At minima/maxima must be zero
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         else
            rDLI=1.d0/rDL
            rRLL=rDLL*rDLI
            rRL =rDM *rDLI
            rBetaL=(-2.d0*rRLL+11.d0+rRL*24.d0-rDR*rDLI*3.d0)*rC130
            rLimL=MAX(0.d0,MIN(2.d0,2.d0*rRL,rBetaL))
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLimL*rDL
         end if

         if(abs(rDR).lt.rSmall)then !.eq.0.d0)then
            !BetaR will be infinite giving rLimL=2
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
         else
            rDRI=1.d0/rDR
            rRR =rDM *rDRI
            rRRR=rDRR*rDRI
            rBetaR=(-2.d0*rRRR+11.d0+rRR*24.d0-rDL*rDRI*3.d0)*rC130
            rLimR=MAX(0.d0,MIN(2.d0,2.d0*rRR,rBetaR))
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLimR*rDR
         end if

      END DO

      if(iLocalRecons.eq.1)then
         v=iLocalNvar-1

         rDLL=rStencil(v,2)-rStencil(v,1) !LL-LLL - i-3/2
         rDL=rStencil(v,3)-rStencil(v,2)  !L-LL - i-1/2
         rDM=rStencil(v,4)-rStencil(v,3)  !R-L - i+1/2
         rDR=rStencil(v,5)-rStencil(v,4)  !RR-R - i+3/2
         rDRR=rStencil(v,6)-rStencil(v,5) !RRR-RR - i+5/2

         if(abs(rDL).lt.rSmall)then !.eq.0.d0)then
            !BetaL will be infinite giving rLimL=2 OR 0.d0 At minima/maxima must be zero
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)
         else
            rDLI=1.d0/rDL
            rRLL=rDLL*rDLI
            rRL =rDM *rDLI
            rBetaL=(-2.d0*rRLL+11.d0+rRL*24.d0-rDR*rDLI*3.d0)*rC130
            rLimL=MAX(0.d0,MIN(2.d0,2.d0*rRL,rBetaL))
            rLimitedVariables(v,1)=rStencil(v,iHalfStencilSize)+0.5d0*rLimL*rDL
         end if

         if(abs(rDR).lt.rSmall)then !.eq.0.d0)then
            !BetaR will be infinite giving rLimL=2
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)
         else
            rDRI=1.d0/rDR
            rRR =rDM *rDRI
            rRRR=rDRR*rDRI
            rBetaR=(-2.d0*rRRR+11.d0+rRR*24.d0-rDL*rDRI*3.d0)*rC130
            rLimR=MAX(0.d0,MIN(2.d0,2.d0*rRR,rBetaR))
            rLimitedVariables(v,2)=rStencil(v,iHalfStencilSize+1)-0.5d0*rLimR*rDR
         end if
      end if

end select

end Subroutine VariableReconstruction

!> subroutine responsible for computing cell interface primitive and conserved variables for use in the Riemann solver (left and right quantities). 
!> The parts of the routine used depend on whether the initial reconstruction was for primitive or conservative variables. This populates total energy if the primitive variables are reconstructed, and computes pressures/temperature/gamma if conserved variables are used. 
!> The key input to this subroutine is rLimitedVariables which has either limited conserved variables or limited primitive variables. This is then modified, and used to computed rLimitedPrimitiveVariables=[u,v,w,rho i,rho,g,p,T] and rLimitedVariables=[rho u, rho v, rho w,rho e, rho,g,p,T] for single species.
Subroutine ComputePCVariables(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant,iLocalRecons)
implicit none
!>Local array sizes
integer, intent(in):: iLocalNVar,iLocalNumberofSpecies,iLocalEquationOfState
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(inout)::rLimitedVariables
!> Left and right limited primitive variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(out)::rLimitedPrimitiveVariables
!> Left and right limited kinetic energy
real(8),dimension (2),intent(out)::rLimitedKineticEnergy
!>left (=1) or right (=2) interface index
!> Local gas properties
real(8),dimension (iLocalNumberOfSpecies),intent(in)::rLocalSpecificGasConstant
real(8),dimension (iLocalNumberOfSpecies,2,5),intent(in)::rLocalCv,rLocalCp
!> Flag for the Reconstruction approach
integer,intent(in)::iLocalRecons
!> Integer to indicate left/right side of interface
integer:: iSide
!> One over density
real(8):: rInverseDensity
!> 1/2
real(8):: rC12,rC13,rC14,rC15
!>variable index
integer:: v
!>loop variable
integer::iLoop
!>1/(gi-1)
real(8),dimension(iLocalNumberOfSpecies):: rInverseGammaM1
!> Density of each species
real(8),dimension(iLocalNumberOfSpecies)::rLocalRho
!>Sum(zi/(gi-1))
real(8) rVolFracInverseGammaM1
!> Specific gas constant of the mixture
real(8) rSpecificGasConstantMixture
!> Temperature of each species and it's powers
real(8) alpha,alpha2,alpha3,alpha4,alpha5
!> Difference at 1000 k
real(8),dimension(iLocalNumberOfSpecies):: Diff_at_1000
!> One over the specific gas constant
real(8),dimension(iLocalNumberOfSpecies):: rLocalSpecificGasConstantI
!> Gamma of each species
real(8) Gammai
!> Intermediate sum
real(8) sumGammai
!> Van der Waals constants of the mixture
real(8) rVDWaMixture,rVDWbMixture
!> Coefficients in the Van Der Waals case
real(8),dimension(iLocalNumberOfSpecies) :: rLocalCa,rLocalCb
!> 1/(specific gas constant)
real(8) rInverseSpecificGasConstant
!> Auxiliary sum to calculate density of the mixture
real(8) rSumRho
!> The integral c_v from 0 to 298.15K and RT_0/W_k (T_{0,k}=\rho_k \_k (int_0^T_0 Cv dT + RT_0/W_k) in Flamenco guide)
real(8),dimension(iLocalNumberOfSpecies):: rT0k
!> T_0 - reference temperature and it's powers for definition of internal energy
real(8):: rT0,rT02,rT03,rT04,rT05
!> Cp of the mixture
real(8) rCpMixture
!> Cv of the mixture
real(8) rCvMixture

select case(iLocalNumberOfSpecies)

   case(1)

      !Single species
      select case (iLocalEquationOfState)
            !> Single species Ideal gas Fixed gamma (perfect gas)
         case(1)
            if (iLocalRecons.eq.0) then
               ! Compute primitive variables from the conserved quantities
               rC12=0.5d0
               rInverseSpecificGasConstant=1.d0/rLocalSpecificGasConstant(1)
               rInverseGammaM1(1)= rLocalCv(1,1,1)*rInverseSpecificGasConstant
               !Gamma
               rLimitedPrimitiveVariables(iLocalNVar-2,:) = 1.d0+1.d0/rInverseGammaM1(1)
               do iSide=1,2
                  ! Velocities
                  rInverseDensity=1.d0/rLimitedVariables(iLocalNVar-3,iSide)
                  do v=1,3
                     rLimitedPrimitiveVariables(v,iSide)=rLimitedVariables(v,iSide)*rInverseDensity
                  end do
                  !Calculate KE, IE, Pressure, Gamma and Temperature
                  rLimitedKineticEnergy(iSide)=rC12*rLimitedVariables(iLocalNVar-3,iSide)*(rLimitedPrimitiveVariables(1,iSide)**2+rLimitedPrimitiveVariables(2,iSide)**2+rLimitedPrimitiveVariables(3,iSide)**2)
                  !>Internal energy
                  rLimitedPrimitiveVariables(4,iSide)=rLimitedVariables(4,iSide)-rLimitedKineticEnergy(iSide)
                  !>Pressure
                  rLimitedPrimitiveVariables(iLocalNVar-1,iSide)=rLimitedPrimitiveVariables(4,iSide)*(rLimitedPrimitiveVariables(iLocalNVar-2,iSide)-1.d0)
                  !Compute temperature from pressure
                  rLimitedPrimitiveVariables(iLocalNVar,iSide) = rInverseDensity*rLimitedPrimitiveVariables(iLocalNVar-1,iSide)*rInverseSpecificGasConstant
               end do
               !Density for primitive variable array
               rLimitedPrimitiveVariables(iLocalNVar-3,:)=rLimitedVariables(iLocalNVar-3,:)
               !Copy computed density, pressure, gamma and temperature back to rLimitedVariables for completeness
               do iLoop=iLocalNVar-2,iLocalNVar
                  rLimitedVariables(iLoop,:)=rLimitedPrimitiveVariables(iLoop,:)
               end do
            else
               !compute conserved variables from primitive
               rC12=0.5d0
               rInverseSpecificGasConstant=1.d0/rLocalSpecificGasConstant(1)
               rInverseGammaM1(1)= rLocalCv(1,1,1)*rInverseSpecificGasConstant
               !Gamma
               rLimitedPrimitiveVariables(iLocalNVar-2,:) = 1.d0+1.d0/rInverseGammaM1(1)
               !Put limited pressure into the Primitive variable array
               rLimitedPrimitiveVariables(iLocalNVar-1,:) =rLimitedVariables(iLocalNVar-1,:)
               !Density for primitive variable array
               rLimitedPrimitiveVariables(iLocalNVar-3,:)=rLimitedVariables(iLocalNVar-3,:)
               do iSide=1,2
                  ! Compute conserved variables from the primitive quantities
                  do v=1,3
                     !first put limited velocities in the primitive variable array
                     rLimitedPrimitiveVariables(v,iSide)=rLimitedVariables(v,iSide)
                     !momenta in conserved array
                     rLimitedVariables(v,iSide)=rLimitedPrimitiveVariables(v,iSide)*rLimitedVariables(iLocalNVar-3,iside)
                  end do
                  rInverseDensity=1.d0/rLimitedVariables(iLocalNVar-3,iSide)
                  !Calculate KE, IE, Temperature
                  rLimitedKineticEnergy(iSide)=rC12*rLimitedVariables(iLocalNVar-3,iSide)*(rLimitedPrimitiveVariables(1,iSide)**2+rLimitedPrimitiveVariables(2,iSide)**2+rLimitedPrimitiveVariables(3,iSide)**2)
                  !>Internal energy
                  rLimitedPrimitiveVariables(4,iSide)=rLimitedPrimitiveVariables(iLocalNVar-1,iSide)*rInverseGammaM1(1)
                  !Compute temperature from pressure
                  rLimitedPrimitiveVariables(iLocalNVar,iSide) = rInverseDensity*rLimitedPrimitiveVariables(iLocalNVar-1,iSide)*rInverseSpecificGasConstant
                  !>total energy
                  rLimitedVariables(4,iSide)=rLimitedKineticEnergy(iSide)+rLimitedPrimitiveVariables(4,iSide)
               end do
               !Copy computed gamma and temperature back to rLimitedVariables for completeness
               rLimitedVariables(iLocalNVar-2,:)=rLimitedPrimitiveVariables(iLocalNVar-2,:)
               rLimitedVariables(iLocalNVar  ,:)=rLimitedPrimitiveVariables(iLocalNVar,:)
            end if
         case default
            ! call Exception("Single gas primitive var calc only implemented for fixed gamma ideal gases", 1)
            print*,'Single gas primitive var calc only implemented for fixed gamma ideal gases'
            stop
      end select

   case default

      rC12=0.5d0
      rC13=1.d0/3.d0
      rC14=0.25d0
      rC15=0.2d0

      if (iLocalRecons.eq.0) then
         do iSide=1,2
            ! Compute primitive variables from the conserved quantities
            ! Density
            rSumRho = 0.d0
            do iLoop=1,iLocalNumberOfSpecies
               rSumRho = rSumRho + rLimitedVariables(4+iLoop,iSide)
            end do
            rLimitedPrimitiveVariables(iLocalNVar-3,iSide) = rSumRho
            rLimitedVariables(iLocalNVar-3,iSide) = rSumRho
            ! Velocities
            rInverseDensity=1.d0/rSumRho
            do v=1,3
               rLimitedPrimitiveVariables(v,iSide)=rLimitedVariables(v,iSide)*rInverseDensity
            end do
            ! Compute Y_i
            do v=5,4+iLocalNumberOfSpecies
               rLimitedPrimitiveVariables(v,iSide)=rLimitedVariables(v,iSide)*rInverseDensity
            end do
            !Calculate KE, IE, Pressure, Gamma and Temperature
            rLimitedKineticEnergy(iSide)=rC12*rLimitedVariables(iLocalNVar-3,iSide)*(rLimitedPrimitiveVariables(1,iSide)**2+rLimitedPrimitiveVariables(2,iSide)**2+rLimitedPrimitiveVariables(3,iSide)**2)
            rLimitedPrimitiveVariables(4,iSide)=rLimitedVariables(4,iSide)-rLimitedKineticEnergy(iSide)
            !Compute the remaining variables through internal energy
            call ComputePrimitiveFromReconstructed(iSide,rLimitedVariables,rLimitedPrimitiveVariables,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant)
         end do
         !Copy computed pressure, gamma and temperature back to rLimitedVariables for completeness
         do iLoop=iLocalNVar-2,iLocalNVar
            rLimitedVariables(iLoop,:)=rLimitedPrimitiveVariables(iLoop,:)
         end do
      else
         do iSide=1,2
            ! Compute conserved variables from primitive
            !Recover rho_i z_i  and compute density
            rSumRho = 0.d0
            do iLoop=1,iLocalNumberOfSpecies
               rLimitedPrimitiveVariables(4+iLoop,iSide)=rLimitedVariables(4+iLoop,iSide)
                  rSumRho = rSumRho + rLimitedVariables(4+iLoop,iSide)
            end do
            rLimitedPrimitiveVariables(iLocalNVar-3,iSide) = rSumRho
            rLimitedVariables(iLocalNVar-3,iSide) = rSumRho
            !momenta
            do v=1,3
               rLimitedPrimitiveVariables(v,iSide)=rLimitedVariables(v,iSide)
               rLimitedVariables(v,iSide)=rLimitedVariables(v,iSide)*rLimitedPrimitiveVariables(iLocalNVar-3,iSide)
            end do
            !Put limited pressure into the Primitive variable array
            rLimitedPrimitiveVariables(iLocalNVar-1,iSide) =rLimitedVariables(iLocalNVar-1,iSide)
            !Calculate KE
            rLimitedKineticEnergy(iSide)=rC12*rLimitedVariables(iLocalNVar-3,iSide)*(rLimitedPrimitiveVariables(1,iSide)**2+rLimitedPrimitiveVariables(2,iSide)**2+rLimitedPrimitiveVariables(3,iSide)**2)
            !Compute internal energy, gamma and temperature of the mixture from!> Single species Ideal gas Fixed gamma (perfect gas) reconstructed pressure and partial densities.
            select case (iLocalEquationOfState)
               case(1) !Perfect gas
                  ! Mixture specific heats used to calculate gamma.
                  rCpMixture = 0.d0
                  rCvMixture = 0.d0
                  do iLoop=1,iLocalNumberOfSpecies
                     rCpMixture = rCpMixture + rLocalCp(iLoop,1,1)*rLimitedPrimitiveVariables(4+iLoop,iSide)
                     rCvMixture = rCvMixture + rLocalCv(iLoop,1,1)*rLimitedPrimitiveVariables(4+iLoop,iSide)
                  end do
                  !Gamma
                  rLimitedPrimitiveVariables(iLocalNVar-2,iSide) = rCpMixture/rCvMixture
                  !Internal energy
                  rLimitedPrimitiveVariables(4,iSide)=rLimitedPrimitiveVariables(iLocalNVar-1,iSide)/(rLimitedPrimitiveVariables(iLocalNVar-2,iSide)-1.d0)
                  !>total energy
                  rLimitedVariables(4,iSide)=rLimitedKineticEnergy(iSide)+rLimitedPrimitiveVariables(4,iSide)         
                  rInverseDensity = 1.d0/rLimitedPrimitiveVariables(iLocalNVar-3,iSide)
                  !Compute R of the mixture
                  rSpecificGasConstantMixture = 0.d0
                  do iLoop=1,iLocalNumberOfSpecies
                     rSpecificGasConstantMixture = rSpecificGasConstantMixture + rLimitedVariables(4+iLoop,iSide)*rInverseDensity*rLocalSpecificGasConstant(iLoop)
                  end do
                  !Compute temperature from pressure
                  rLimitedPrimitiveVariables(iLocalNVar,iSide) = rInverseDensity*rLimitedPrimitiveVariables(iLocalNVar-1,iSide)/rSpecificGasConstantMixture
               case default
                  ! call Exception("Multispecies primitive var calc only implemented for fixed gamma ideal gases", 1)
                  print*,'Multispecies primitive var calc only implemented for fixed gamma ideal gases'
                  stop

            end select
         end do
         !Copy computed gamma and temperature back to rLimitedVariables for completeness
         rLimitedVariables(iLocalNVar-2,:)=rLimitedPrimitiveVariables(iLocalNVar-2,:)
         rLimitedVariables(iLocalNVar  ,:)=rLimitedPrimitiveVariables(iLocalNVar,:)
      end if
end select

end Subroutine ComputePCVariables

Subroutine ComputePrimitiveFromReconstructed(iSide,rLimitedVariables,rLimitedPrimitiveVariables,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalCv,rLocalCp,rLocalSpecificGasConstant)
implicit none
!>Local array sizes
integer, intent(in):: iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(in)::rLimitedVariables
!> Left and right limited primitive variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(inout)::rLimitedPrimitiveVariables
!> Local gas properties
real(8),dimension (iLocalNumberOfSpecies),intent(in)::rLocalSpecificGasConstant
real(8),dimension (iLocalNumberOfSpecies,2,5),intent(in)::rLocalCv,rLocalCp
!> Left or Right
integer,intent(in)::iSide
!> One over density
real(8) rInverseDensity
!>1/(gi-1)
real(8),dimension(iLocalNumberOfSpecies):: rInverseGammaM1
!> Loop index
integer iLoop
!>Sum(zi/(gi-1))
real(8) rVolFracInverseGammaM1
!> Specific gas constant of the mixture
real(8) rSpecificGasConstantMixture
!> Van der Waals constants of the mixture
real(8) rVDWaMixture,rVDWbMixture
!> Coefficients in the Van Der Waals case
real(8),dimension(iLocalNumberOfSpecies) :: rLocalCa,rLocalCb
!> Temperature of each species
real(8) alpha,alpha2,alpha3,alpha4
!> Gamma of each species
real(8) Gammai
!> Intermediate sum
real(8) sumGammai
!> Density of each species
real(8),dimension(iLocalNumberOfSpecies)::rLocalRho
!>Auxiliar vector for the Newton-Raphson solver
real(8),dimension(iLocalNVar)::rThermodynamicVector
!> Cp of the mixture
real(8) rCpMixture
!> Cv of the mixture
real(8) rCvMixture

select case (iLocalNumberOfSpecies)
   case(1)

      ! call Exception("ComputePrimitiveFromReconstructed should not be called for single gas primitive var calc in inviscid_fluxes.f90", 1)
      print*,'ComputePrimitiveFromReconstructed should not be called for single gas primitive var calc in inviscid_fluxes.f90'
      stop

   case default
      select case (iLocalEquationOfState)
         !> Multispecies. Ideal gas Fixed gamma (perfect gas)
         case(1)

            rInverseDensity=1.d0/rLimitedVariables(iLocalNVar-3,iSide)
            ! Mixutre specific heats used to calculate gamma. Note - actually rho*Cp and rho*Cv
            rCpMixture = 0.d0
            rCvMixture = 0.d0
            do iLoop=1,iLocalNumberOfSpecies
               rCpMixture = rCpMixture + rLocalCp(iLoop,1,1)*rLimitedVariables(4+iLoop,iSide)
               rCvMixture = rCvMixture + rLocalCv(iLoop,1,1)*rLimitedVariables(4+iLoop,iSide)
            end do
            !Gamma
            rLimitedPrimitiveVariables(iLocalNVar-2,iSide) = rCpMixture/rCvMixture
            !Pressure
            rLimitedPrimitiveVariables(iLocalNVar-1,iSide)=(rLimitedPrimitiveVariables(iLocalNVar-2,iSide)-1.d0)*rLimitedPrimitiveVariables(4,iSide)
            !Temperature
            !Compute R of the mixture
            rSpecificGasConstantMixture = 0.d0
            do iLoop=1,iLocalNumberOfSpecies
               rSpecificGasConstantMixture = rSpecificGasConstantMixture + rLimitedVariables(4+iLoop,iSide)*rInverseDensity*rLocalSpecificGasConstant(iLoop)
            end do
            !Compute temperature from pressure
            rLimitedPrimitiveVariables(iLocalNVar,iSide) = rInverseDensity*rLimitedPrimitiveVariables(iLocalNVar-1,iSide)/rSpecificGasConstantMixture
         case default

            ! call Exception("ComputePrimitiveFromReconstructed only implemented for ideal gas EoS", 1)
            print*,'ComputePrimitiveFromReconstructed only implemented for ideal gas EoS'
            stop
      end select
end select

end subroutine ComputePrimitiveFromReconstructed

!>subroutine responsible for computing the low Mach correction. Not used for 1st order in space.
Subroutine LowMachCorrection(rLimitedVariables,rLimitedPrimitiveVariables,rLimitedKineticEnergy,iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState,rLocalSpecificGasConstant,rLocalCv,rLocalCp)
implicit none
!>Local array sizes
integer, intent(in):: iLocalNVar,iLocalNumberOfSpecies,iLocalEquationOfState
!> Local gas properties
real(8),dimension (iLocalNumberOfSpecies),intent(in)::rLocalSpecificGasConstant
real(8),dimension (iLocalNumberOfSpecies,2,5),intent(in)::rLocalCv,rLocalCp
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(inout)::rLimitedVariables
!> Left and right limited primitive variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(inout)::rLimitedPrimitiveVariables
!> Left and right limited kinetic energy
real(8),dimension(2),intent(inout)::rLimitedKineticEnergy
!>variable index
integer v, iLoop
!>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR LEFT SIDE OF RP
real(8)::DL,UL,VL,WL,PL,GL
!>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR RIGHT SIDE OF RP
real(8)::DR,UR,VR,WR,PR,GR
!>specific kinetic energy*2 on the left and right
real(8) KEL,KER
!> left and right Mach number
real(8) ML,MR
!> Low Mach limiter function
real(8) Lim
!>VELOCITY DIFFERENCE
real(8)::DU,DV,DW
!>VELOCITY SUMS
real(8)::SU,SV,SW
!>1/2
real(8) rC12
!>Speed of sound left and right, and squared
real(8):: SOSL,SOSR,SOSL2,SOSR2

rC12=1.d0/2.d0

UL = rLimitedPrimitiveVariables(1,1)
VL = rLimitedPrimitiveVariables(2,1)
WL = rLimitedPrimitiveVariables(3,1)
PL = rLimitedPrimitiveVariables(iLocalNVar-1,1)
DL = rLimitedPrimitiveVariables(iLocalNVar-3,1)
GL = rLimitedPrimitiveVariables(iLocalNVar-2,1)

UR = rLimitedPrimitiveVariables(1,2)
VR = rLimitedPrimitiveVariables(2,2)
WR = rLimitedPrimitiveVariables(3,2)
PR = rLimitedPrimitiveVariables(iLocalNVar-1,2)
DR = rLimitedPrimitiveVariables(iLocalNVar-3,2)
GR = rLimitedPrimitiveVariables(iLocalNVar-2,2)

!Speed of sound from mixture gammas
SOSL2=GL*PL/DL
SOSR2=GR*PR/DR

KEL=UL**2+VL**2+WL**2
ML=SQRT(KEL/SOSL2)
KER=UR**2+VR**2+WR**2
MR=SQRT(KER/SOSR2)

LIM=MIN(MAX(ML,MR),1.d0)

DU=LIM*(UL-UR)
DV=LIM*(VL-VR)
DW=LIM*(WL-WR)

SU=(UL+UR)
SV=(VL+VR)
SW=(WL+WR)

UL=rC12*(SU+DU)
VL=rC12*(SV+DV)
WL=rC12*(SW+DW)
UR=rC12*(SU-DU)
VR=rC12*(SV-DV)
WR=rC12*(SW-DW)

rLimitedPrimitiveVariables(1,1)=UL
rLimitedPrimitiveVariables(2,1)=VL
rLimitedPrimitiveVariables(3,1)=WL
rLimitedPrimitiveVariables(1,2)=UR
rLimitedPrimitiveVariables(2,2)=VR
rLimitedPrimitiveVariables(3,2)=WR

do v=1,3
   rLimitedVariables(v,1)=rLimitedPrimitiveVariables(v,1)*DL
   rLimitedVariables(v,2)=rLimitedPrimitiveVariables(v,2)*DR
end do

rLimitedKineticEnergy(1)=rC12*DL*(UL**2+VL**2+WL**2)
rLimitedKineticEnergy(2)=rC12*DR*(UR**2+VR**2+WR**2)
rLimitedVariables(4,1)=rLimitedPrimitiveVariables(4,1)+rLimitedKineticEnergy(1)
rLimitedVariables(4,2)=rLimitedPrimitiveVariables(4,2)+rLimitedKineticEnergy(2)

end Subroutine LowMachCorrection

!>subroutine responsible for solving the curvilinear Riemann problems. This subroutine utilises the HLLC solver.
Subroutine RiemannSolver(rLimitedVariables,rLimitedPrimitiveVariables,rStencil,iLocalNVar,iLocalNVarComputed,iLocalNcHalo,iLocalNumberOfSpecies,iLocalEquationOfState,iHalfStencilSize,rLocalSpecificGasConstant,rLocalCv,rLocalCp,rNJD1DX,rNJD1DY,rNJD1DZ,rBeta,rLocalFlux,i,j,k,iDir)
implicit none
!>Local array sizes
integer, intent(in):: iLocalNVar,iLocalNVarComputed,iLocalNumberOfSpecies,iLocalEquationOfState,iLocalNcHalo
!> half the stencil size (note that stencils must be an even number of points)
integer iHalfStencilSize
!> Local gas properties
real(8),dimension (iLocalNumberOfSpecies),intent(in)::rLocalSpecificGasConstant
real(8),dimension (iLocalNumberOfSpecies,2,5),intent(in)::rLocalCv,rLocalCp
!> Fluxes for each equation
real(8),dimension (4+iLocalNumberOfSpecies),intent(out)::rLocalFlux
!> Left and right limited variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(in)::rLimitedVariables
!> Left and right limited primitive variables (left=(:,1), right=(:,2))
real(8),dimension (iLocalNVar,2),intent(in)::rLimitedPrimitiveVariables
!> Stencil
real(8),dimension (iLocalNVar,2*iLocalNcHalo)::rStencil
!>THE DERIVATIVES OF XI, ETA, ZETA IN DIRECTION 1 (THE DIRECTION ALONG THE CURRENT LINE OF CELLS) MULTLIPLIED BY THE JACOBIAN J AND NORMALISED TO BECOME UNIT NORMAL VECTOR COMPONENTS
real(8),intent(in)::rNJD1DX,rNJD1DY,rNJD1DZ
!>Metric vector magnitude
real(8),intent(in):: rBeta
!>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR LEFT SIDE OF RP
real(8)::DL,UL,VL,WL,PL,GL,CL
!>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR RIGHT SIDE OF RP
real(8)::DR,UR,VR,WR,PR,GR,CR
!>UL AND UR NORMAL TO THE INTERFACE
real(8)::ULN,URN
!>WAVE SPEEDS FOR THE L, STAR AND RIGHT STATES
real(8):: SL, SM, SR
!>1/(SL-SM)
real(8)::C1SLSM
!>PRESSURE IN THE STAR REGION
real(8)::PSTAR
!>SL-ULN WHERE ULN IS THE FACE NORMAL VELOCITY FROM THE LEFT
real(8)::SLULN
!>PSTAR-PL
real(8)::PSPL
!>SM*C1SLSM
real(8)::SMC1SLSM
!>1/(SR-SM)
real(8)::C1SRSM
!>SR-URN WHERE URN IS THE FACE NORMAL VELOCITY FROM THE RIGHT
real(8)::SRURN
!>PSTAR-PR
real(8):: PSPR
!>SM*C1SRSM
real(8)::SMC1SRSM
!>GAMMA AND RELATED FUNCTIONS
real(8):: G, G1, G2, G3, G4, G5, G6, G7,r2GI,rGM1I,rGP1I
integer i,j,k
!>Loop
integer::iLoop
!> 1.d0/2.d0
real(8):: rC12
!> Direction
integer::iDir


rLocalFlux=0.d0

UL = rLimitedPrimitiveVariables(1,1)
VL = rLimitedPrimitiveVariables(2,1)
WL = rLimitedPrimitiveVariables(3,1)
PL = rLimitedPrimitiveVariables(iLocalNVar-1,1)
GL = rLimitedPrimitiveVariables(iLocalNVar-2,1)
DL = rLimitedPrimitiveVariables(iLocalNVar-3,1)

UR = rLimitedPrimitiveVariables(1,2)
VR = rLimitedPrimitiveVariables(2,2)
WR = rLimitedPrimitiveVariables(3,2)
PR = rLimitedPrimitiveVariables(iLocalNVar-1,2)
GR = rLimitedPrimitiveVariables(iLocalNVar-2,2)
DR = rLimitedPrimitiveVariables(iLocalNVar-3,2)

!> Moving mesh velocities
UL = UL
VL = VL
WL = WL
UR = UR
VR = VR
WR = WR

!Speed of sound from mixture gammas
CL=SQRT(GL*PL/DL)
CR=SQRT(GR*PR/DR)

!CELL INTERFACE NORMAL VELOCITIES
ULN = (UL*rNJD1DX+VL*rNJD1DY+WL*rNJD1DZ)
URN = (UR*rNJD1DX+VR*rNJD1DY+WR*rNJD1DZ)

!C        Compute fluxes FDL and FDR at CDL and CDR
!C        Calculate estimates for wave speeds using adaptive
!C        approximatd-state Riemann solvers
!C        Use average gamma for wavespeed calculations
rC12=1.d0/2.d0
G  = (GL+GR)*rC12
r2GI=rC12/G
rGM1I=1.d0/(G - 1.0d0)
rGP1I=1.d0/(G + 1.0d0)
G1 = (G - 1.0d0)*r2GI
G2 = (G + 1.0d0)*r2GI
G3 = 2.0d0*G*rGM1I
G4 = 2.0d0*rGM1I
G5 = 2.0d0*rGP1I
G6 = (G - 1.0d0)*rGP1I
G7 = (G - 1.0d0)*rC12

!Estimate the wavespeeds
CALL ESTIME(SL,SM,SR,G1,G2,G3,G4,G5,G6,G7,DL,PL,CL,DR,PR,CR,ULN,URN)

IF(SL.GE.0.0d0)THEN

   !C           Right-going supersonic flow

   rLocalFlux(1) = rLimitedVariables(1,1)*ULN + PL*rNJD1DX
   rLocalFlux(2) = rLimitedVariables(2,1)*ULN + PL*rNJD1DY
   rLocalFlux(3) = rLimitedVariables(3,1)*ULN + PL*rNJD1DZ
   rLocalFlux(4) = ULN*(rLimitedVariables(4,1) + PL)
   if(iLocalNumberOfSpecies.eq.1)then
      rLocalFlux(iLocalNVar-3) = DL*ULN
   else
      do iLoop=1,iLocalNumberOfSpecies
         rLocalFlux(4+iLoop) = rLimitedVariables(4+iLoop,1)*ULN
      end do
   end if

   !REMULTIPLY BY THE MAGNITUDE OF THE VECTOR TO RECOVER ORIGINAL JD1DZ'S, AND THE SIGN OF THE JACOBIAN
   rLocalFlux(1:4+iLocalNumberOfSpecies) = rLocalFlux(1:4+iLocalNumberOfSpecies)*rBeta

ENDIF

IF(SL.LE.0.0d0.AND.SR.GE.0.0d0)THEN

   !!$        !C           Subsonic flow

   IF(SM.GE.0.0d0)THEN

      !C              Subsonic flow to the right

      C1SLSM=1.d0/(SL - SM)
      PSTAR=DL*(ULN-SL)*(ULN-SM)+PL
      SLULN=SL-ULN
      PSPL=PSTAR-PL
      SMC1SLSM=SM*C1SLSM

      rLocalFlux(1) = (rLimitedVariables(1,1)*SLULN+PSPL*rNJD1DX)*SMC1SLSM+PSTAR*rNJD1DX
      rLocalFlux(2) = (rLimitedVariables(2,1)*SLULN+PSPL*rNJD1DY)*SMC1SLSM+PSTAR*rNJD1DY
      rLocalFlux(3) = (rLimitedVariables(3,1)*SLULN+PSPL*rNJD1DZ)*SMC1SLSM+PSTAR*rNJD1DZ
      rLocalFlux(4) = (rLimitedVariables(4,1)*SLULN+PSTAR*SL-PL*ULN)*SMC1SLSM
      if(iLocalNumberOfSpecies.eq.1)then
         rLocalFlux(iLocalNVar-3) = DL*SLULN*SMC1SLSM
      else
         do iLoop=1,iLocalNumberOfSpecies
            rLocalFlux(4+iLoop) = rLimitedVariables(4+iLoop,1)*SLULN*SMC1SLSM
         end do
      end if

      !REMULTIPLY BY THE MAGNITUDE OF THE VECTOR TO RECOVER ORIGINAL JD1DZ'S, AND THE SIGN OF THE JACOBIAN
      rLocalFlux(1:4+iLocalNumberOfSpecies)=rLocalFlux(1:4+iLocalNumberOfSpecies)*rBeta

   ELSE

      !C              Subsonic flow to the left

      C1SRSM=1.d0/(SR - SM)
      PSTAR=DR*(URN-SR)*(URN-SM)+PR
      SRURN=SR-URN
      PSPR=PSTAR-PR
      SMC1SRSM=SM*C1SRSM

      rLocalFlux(1) = (rLimitedVariables(1,2)*SRURN+PSPR*rNJD1DX)*SMC1SRSM+PSTAR*rNJD1DX
      rLocalFlux(2) = (rLimitedVariables(2,2)*SRURN+PSPR*rNJD1DY)*SMC1SRSM+PSTAR*rNJD1DY
      rLocalFlux(3) = (rLimitedVariables(3,2)*SRURN+PSPR*rNJD1DZ)*SMC1SRSM+PSTAR*rNJD1DZ
      rLocalFlux(4) = (rLimitedVariables(4,2)*SRURN+PSTAR*SR-PR*URN)*SMC1SRSM
      if(iLocalNumberOfSpecies.eq.1)then
         rLocalFlux(iLocalNVar-3) = DR*SRURN*SMC1SRSM
      else
         do iLoop=1,iLocalNumberOfSpecies
            rLocalFlux(4+iLoop) = rLimitedVariables(4+iLoop,2)*SRURN*SMC1SRSM
         end do
      end if

      !REMULTIPLY BY THE MAGNITUDE OF THE VECTOR TO RECOVER ORIGINAL JD1DZ'S, AND THE SIGN OF THE JACOBIAN
      rLocalFlux(1:4+iLocalNumberOfSpecies)=rLocalFlux(1:4+iLocalNumberOfSpecies)*rBeta

   ENDIF
ENDIF

IF(SR.LE.0.0d0)THEN

   !C           Left-going supersonic flow

   rLocalFlux(1) = rLimitedVariables(1,2)*URN + PR*rNJD1DX
   rLocalFlux(2) = rLimitedVariables(2,2)*URN + PR*rNJD1DY
   rLocalFlux(3) = rLimitedVariables(3,2)*URN + PR*rNJD1DZ
   rLocalFlux(4) = URN*(rLimitedVariables(4,2) + PR)
   if(iLocalNumberOfSpecies.eq.1)then
      rLocalFlux(iLocalNVar-3) = DR*URN
   else
      do iLoop=1,iLocalNumberOfSpecies
         rLocalFlux(4+iLoop) = rLimitedVariables(4+iLoop,2)*URN
      end do
   end if

   !REMULTIPLY BY THE MAGNITUDE OF THE VECTOR TO RECOVER ORIGINAL JD1DZ'S, AND THE SIGN OF THE JACOBIAN
   rLocalFlux(1:4+iLocalNumberOfSpecies) = rLocalFlux(1:4+iLocalNumberOfSpecies)*rBeta

ENDIF

end Subroutine RiemannSolver
!>     Purpose: to compute wave speed estimates for the HLLC Riemann
!>              solver using and adaptive approximate-state Riemann
!>              solver including the PVRS, TRRS and TSRS solvers.
!>              See Chap. 9, Ref. 1 Toro.
SUBROUTINE ESTIME(SL,SM,SR,G1,G2,G3,G4,G5,G6,G7,DL,PL,CL,DR,PR,CR,ULN,URN)

  !     Declaration of variables

  !>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR LEFT SIDE OF RP
real(8),intent(in)::DL,PL,CL
!>DENSITY,VELS,PRESSURE, GAMMA, SOS FOR RIGHT SIDE OF RP
real(8),intent(in)::DR,PR,CR
!>Inverse sound speeds
real(8) ::CLI,CRI
!>Signal speeds
real(8),intent(out)::SL, SM, SR
!>variables used only within this routine
real(8)    CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,PTL, PTR, QMAX, QUSER, UM
!>GAMMA AND RELATED FUNCTIONS
real(8):: G1, G2, G3, G4, G5, G6, G7
!>UL AND UR NORMAL TO THE INTERFACE
real(8)::ULN,URN

QUSER = 2.0d0

!C     Compute guess pressure from PVRS Riemann solver

CUP  = 0.25d0*(DL + DR)*(CL + CR)
PPV  = 0.5d0*(PL + PR) + 0.5d0*(ULN - URN)*CUP
PPV  = MAX(0.0d0, PPV)
PMIN = MIN(PL,  PR)
PMAX = MAX(PL,  PR)
QMAX = PMAX/PMIN

IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN

    !C        Select PRVS Riemann solver

   PM = PPV
   UM = 0.5d0*(ULN + URN) + 0.5d0*(PL - PR)/CUP

ELSE
   IF(PPV.LT.PMIN)THEN

      !C           Select Two-Rarefaction Riemann solver

      PQ  = (PL/PR)**G1
      CLI=1.d0/CL
      CRI=1.d0/CR
      UM  = (PQ*ULN*CLI + URN*CRI + G4*(PQ - 1.0d0))/(PQ*CLI + 1.0d0*CRI)
      PTL = 1.0d0 + G7*(ULN - UM)*CLI
      PTR = 1.0d0 + G7*(UM - URN)*CRI
      PM  = 0.5d0*(PL*PTL**G3 + PR*PTR**G3)

   ELSE

      !C           Use Two-Shock Riemann solver with PVRS as estimate

      GEL = SQRT(G5/(DL*(G6*PL + PPV)))
      GER = SQRT(G5/(DR*(G6*PR + PPV)))
      PM  = (GEL*PL + GER*PR - (URN - ULN))/(GEL + GER)
      UM  = 0.5d0*(ULN + URN) + 0.5d0*(GER*(PM - PR) - GEL*(PM - PL))
   ENDIF
ENDIF

!C     Find speeds

IF(PM.LE.PL)THEN
   SL = ULN - CL
ELSE
   SL = ULN - CL*SQRT(1.d00 + G2*(PM/PL - 1.0d0))
ENDIF

SM = UM

IF(PM.LE.PR)THEN
   SR = URN + CR
ELSE
   SR = URN + CR*SQRT(1.d00 + G2*(PM/PR - 1.0d0))
ENDIF

END SUBROUTINE ESTIME

subroutine Construct_Cartesian_Xi_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
   implicit none
   !> variable indices
   integer, intent(in) :: i,j,k
   integer j1,k1
   !>Local array sizes
   integer, intent(in):: iLocalNnHalo,iLocalNnx,iLocalNny,iLocalNnz
   !>Local X array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalX
   !>Local Y array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalY
   !>Local Z array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalZ
   !>Vector magnitude
   real(8),intent(out):: rBeta
   !>THE DERIVATIVES OF X IN J DIRECTION
   real(8) rDXD2,rDYD2,rDZD2
   !>THE DERIVATIVES OF X IN K DIRECTION
   real(8) rDXD3,rDYD3,rDZD3
   !>1/2
   real(8) rC12
   !>co-ordinates of the four points
   real(8) rX1,rX2,rX3,rX4,rY1,rY2,rY3,rY4,rZ1,rZ2,rZ3,rZ4
 
   rC12=1.d0/2.d0
 
   !XI DIRECTION - NEED DXIDX,DXIDY,DXIDZ
   !D1=XI=I
   !D2=ETA=J
   !D3=ZETA=K
 
   !COMPUTE DERIVATIVES TO CONSTRUCT THE COMPONENTS OF THE JACOBIAN
 
   J1=J+1
   K1=K+1
 
   rY2=rLocalY(I,J1,K1)
   rZ2=rLocalZ(I,J1,K1)
   rY3=rLocalY(I,J,K)
   rZ3=rLocalZ(I,J,K)
 
   rDYD2=(rY2-rY3)
   rDZD3=(rZ2-rZ3)
 
   rBETA= abs(rDYD2*rDZD3)

end subroutine Construct_Cartesian_Xi_Face_Metrics

subroutine Construct_Cartesian_Et_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
   implicit none
   !> variable indices
   integer, intent(in) :: i,j,k
   integer i1,k1
   !>Local array sizes
   integer, intent(in):: iLocalNnHalo,iLocalNnx,iLocalNny,iLocalNnz
   !>Local X array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalX
   !>Local Y array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalY
   !>Local Z array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalZ
   !>Vector magnitude
   real(8),intent(out):: rBeta
   !>THE DERIVATIVES OF X IN J DIRECTION
   real(8) rDXD2,rDYD2,rDZD2
   !>THE DERIVATIVES OF X IN K DIRECTION
   real(8) rDXD3,rDYD3,rDZD3
   !>1/2
   real(8) rC12
   !>co-ordinates of the four points
   real(8) rX1,rX2,rX3,rX4,rY1,rY2,rY3,rY4,rZ1,rZ2,rZ3,rZ4
 
   rC12=1./2.
 
   !ETA DIRECTION - NEED DETADX,DETADY,DETADZ
   !D1=ETA=J
   !D2=XI=I
   !D3=ZETA=K
   !COMPUTE DERIVATIVES TO CONSTRUCT THE COMPONENTS OF THE JACOBIAN
 
   I1=I+1
   K1=K+1
 
   rX1=rLocalX(I1,J,K)
   rZ1=rLocalZ(I1,J,K)
   rX4=rLocalX(I,J,K1)
   rZ4=rLocalZ(I,J,K1)
 
   rDXD2=(rX1-rX4)
   rDZD3=(rZ4-rZ1)
 
   !COMPUTE D1DX ETC... SECOND ROW OF THE JACOBIAN MATRIX
   !USING CONSERVATIVE DIFFERENCES AS OUTLINED ON P.1034, THOMAS AND LOMBARD, AIAA J. 1979

   rBETA=abs(rDXD2*rDZD3)

end subroutine Construct_Cartesian_Et_Face_Metrics

subroutine Construct_Cartesian_Ze_Face_Metrics(I,J,K,iLocalNnx,iLocalNny,iLocalNnz,iLocalNnHalo,rLocalX,rLocalY,rLocalZ,rBeta)
   implicit none
   !> variable indices
   integer, intent(in) :: i,j,k
   integer i1,j1
   !>Local array sizes
   integer, intent(in):: iLocalNnHalo,iLocalNnx,iLocalNny,iLocalNnz
   !>Local X array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalX
   !>Local Y array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalY
   !>Local Z array
   real(8), dimension (1-iLocalNnHalo:iLocalNnx+iLocalNnHalo,1-iLocalNnHalo:iLocalNny+iLocalNnHalo,1-iLocalNnHalo:iLocalNnz+iLocalNnHalo), intent(in)  :: rLocalZ
   !>Vector magnitude
   real(8),intent(out):: rBeta
   !>THE DERIVATIVES OF X IN J DIRECTION
   real(8) rDXD2,rDYD2,rDZD2
   !>THE DERIVATIVES OF X IN K DIRECTION
   real(8) rDXD3,rDYD3,rDZD3
   !>1/2
   real(8) rC12
   !>co-ordinates of the four points
   real(8) rX1,rX2,rX3,rX4,rY1,rY2,rY3,rY4,rZ1,rZ2,rZ3,rZ4
 
 
   rC12=1./2.
 
   !ETA DIRECTION - NEED DETADX,DETADY,DETADZ
   !D1=ZETA=K
   !D2=XI=I
   !D3=ETA=J
   !COMPUTE DERIVATIVES TO CONSTRUCT THE COMPONENTS OF THE JACOBIAN
 
   I1=I+1
   J1=J+1
 
   rX2=rLocalX(I1,J1,K)
   rY2=rLocalY(I1,J1,K)
   rX3=rLocalX(I,J,K)
   rY3=rLocalY(I,J,K)
 
   rDXD2=(rX2-rX3)
   rDYD3=(rY2-rY3)

   rBETA=abs(rDXD2*rDYD3)

end subroutine Construct_Cartesian_Ze_Face_Metrics
 
                                  

