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
module LES_compact_operators
  USE iso_c_binding
  use LES_objects, only : compact_ops=>compact_ptr ! , mesh_data=>mesh_ptr
  use nvtx
!  use LES_compact, only : compact_ops ! , mesh_data
  IMPLICIT NONE
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double

contains

  subroutine d1x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      dv = zero
      return
    endif
    call nvtxStartRange("d1x")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d1x(2)%ar) ) iop = 2
 !     if( compact_ops%d1x(iop)%id == 0 ) print *,'using compact d1x(',iop,')'
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d1x(iop)%evalx(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dx)
    ! apply metric here or later? here d is scalar or size(dv,1)
    !!$acc parallel loop collapse(2) copy(dv) copyin(compact_ops)
    !do j=1,size(dv,2)
    !do k=1, size(dv,2)
    !   dv(:,j,k) = dv(:,j,k)/compact_ops%dx
    !enddo
    !enddo
    !forall(j=1:size(dv,2),k=1:size(dv,3)) dv(:,j,k) = dv(:,j,k)/compact_ops%dx
    ! dv = dv/mesh_data%d1  ! 3D metric
    call nvtxEndRange()
  end subroutine d1x

  subroutine d1y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      dv = zero
      return
    endif
    call nvtxStartRange("d1y")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d1y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d1y(iop)%evaly(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dy**2)
    ! apply metric here or later? here d is scalar or size(dv,2)
    !forall(i=1:size(dv,1),k=1:size(dv,3)) dv(i,:,k) = dv(i,:,k)/compact_ops%dy
    ! dv = dv/mesh_data%d2  ! 3D metric
    call nvtxEndRange()
  end subroutine d1y

  subroutine d1z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      dv = zero
      return
    endif
    ! symmetry option
    call nvtxStartRange("d1z")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d1z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d1z(iop)%evalz(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dz)
    ! apply metric here or later? here d is scalar or size(dv,3)
    !forall(i=1:size(dv,1),j=1:size(dv,2)) dv(i,j,:) = dv(i,j,:)/compact_ops%dz
    ! dv = dv/mesh_data%d3  ! 3D metric
    call nvtxEndRange()
  end subroutine d1z

  subroutine d2x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    call nvtxStartRange("d2x: dv init")
    if( compact_ops%control%null_opx ) then
      dv = zero
      return
    endif
    call nvtxEndRange
    call nvtxStartRange("d2x")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d2x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call nvtxStartRange("d2x_evalx")
    CALL compact_ops%d2x(iop)%evalx(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dx**2)
    call nvtxEndRange
    ! apply metric here or later? here d is scalar or size(dv,1)
    !!$acc parallel loop collapse(2) copyin(compact_ops) copy(dv)
    !do j=1,size(dv,2)
    !do k=1,size(dv,3)
    !  dv(:,j,k) = dv(:,j,k)/compact_ops%dx**2
    !enddo
    !enddo
    !forall(j=1:size(dv,2),k=1:size(dv,3)) dv(:,j,k) = dv(:,j,k)/compact_ops%dx**2
    ! dv = dv/mesh_data%d1**2  ! 3D metric
    call nvtxEndRange()
  end subroutine d2x

  subroutine d2y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    call nvtxStartRange("d2y: dv init")
    if( compact_ops%control%null_opy ) then
      dv = zero
      return
    endif
    call nvtxEndRange
    call nvtxStartRange("d2y")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d2y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d2y(iop)%evaly(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dy**2)
    ! apply metric here or later? here d is scalar or size(dv,2)
    !!$acc parallel loop collapse(3) copy(dv) copyin(compact_ops)
    !do k=1,size(dv,3)
    !do j=1,size(dv,2)
    !do i=1,size(dv,1)
    !   dv(i,j,k) = dv(i,j,k)/compact_ops%dy**2
    !enddo
    !enddo
    !enddo
    !forall(i=1:size(dv,1),k=1:size(dv,3)) dv(i,:,k) = dv(i,:,k)/compact_ops%dy**2
    ! dv = dv/mesh_data%d2**2  ! 3D metric
    call nvtxEndRange()
  end subroutine d2y

  subroutine d2z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    call nvtxStartRange("d2z: dv init")
    if( compact_ops%control%null_opz ) then
      dv = zero
      return
    endif
    call nvtxEndRange
    call nvtxStartRange("d2z")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d2z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d2z(iop)%evalz(dv,v,vb1,vb2,scalefac_=1.0d0/compact_ops%dz**2)
    ! apply metric here or later? here d is scalar or size(dv,3)
    !!$acc parallel loop collapse(3) copy(dv) copyin(compact_ops)
    !do k=1,size(dv,3)
    !do j=1,size(dv,2)
    !do i=1,size(dv,1)
    !   dv(i,j,k) = dv(i,j,k)/compact_ops%dz**2
    !enddo
    !enddo
    !enddo
    !forall(i=1:size(dv,1),j=1:size(dv,2)) dv(i,j,:) = dv(i,j,:)/compact_ops%dz**2
    ! dv = dv/mesh_data%d3**2  ! 3D metric
    call nvtxEndRange()
  end subroutine d2z

  subroutine d8x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      dv = zero
      return
    endif
    call nvtxStartRange("d8x")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d8x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d8x(iop)%evalx(dv,v,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,1)
 !   forall(j=1:size(dv,2),k=1:size(dv,3)) dv(:,j,k) = dv(:,j,k)/compact_ops%dx**8
    ! dv = dv/mesh_data%d1**8  ! 3D metric
    call nvtxEndRange()
  end subroutine d8x

  subroutine d8y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      dv = zero
      return
    endif
    call nvtxStartRange("d8y")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d8y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d8y(iop)%evaly(dv,v,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,2)
 !   forall(i=1:size(dv,1),k=1:size(dv,3)) dv(i,:,k) = dv(i,:,k)/compact_ops%dy**8
    ! dv = dv/mesh_data%d2**8  ! 3D metric
    call nvtxEndRange()
  end subroutine d8y

  subroutine d8z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      dv = zero
      return
    endif
    call nvtxStartRange("d8z")
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%d8z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    CALL compact_ops%d8z(iop)%evalz(dv,v,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,3)
!    forall(i=1:size(dv,1),j=1:size(dv,2)) dv(i,j,:) = dv(i,j,:)/compact_ops%dz**8
    ! dv = dv/mesh_data%d3**8  ! 3D metric
    call nvtxEndRange()
  end subroutine d8z

  subroutine filterx(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      fv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("filterx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. compact_ops%nop(1) > 1 ) iop = 2
    endif
    ! select filter
    if( nf == compact_ops%control%gfspec ) then
      CALL compact_ops%gfx(iop)%evalx(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%sfspec ) then
      CALL compact_ops%sfx(iop)%evalx(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%tfspec ) then
      CALL compact_ops%tfx(iop)%evalx(fv,v,vb1,vb2)
      return
    endif
    call nvtxEndRange()
  end subroutine filterx

  subroutine filtery(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      fv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("filtery")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. compact_ops%nop(2) > 1 ) iop = 2
    endif
    ! select filter
    if( nf == compact_ops%control%gfspec ) then
      CALL compact_ops%gfy(iop)%evaly(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%sfspec ) then
      CALL compact_ops%sfy(iop)%evaly(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%tfspec ) then
      CALL compact_ops%tfy(iop)%evaly(fv,v,vb1,vb2)
      return
    endif
    call nvtxEndRange()
  end subroutine filtery

  subroutine filterz(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      fv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("filterz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. compact_ops%nop(3) > 1 ) iop = 2
    endif
    ! select filter
    if( nf == compact_ops%control%gfspec ) then
      CALL compact_ops%gfz(iop)%evalz(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%sfspec ) then
      CALL compact_ops%sfz(iop)%evalz(fv,v,vb1,vb2)
      return
    endif
    if( nf == compact_ops%control%tfspec ) then
      CALL compact_ops%tfz(iop)%evalz(fv,v,vb1,vb2)
      return
    endif
    call nvtxEndRange()
  end subroutine filterz

  subroutine islx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opx ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("islx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%islx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%islx(iop)%evalx(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine islx

  subroutine isly(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opy ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("isly")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%isly(2)%ar) ) iop = 2
    endif
    CALL compact_ops%isly(iop)%evaly(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine isly

  subroutine islz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opz ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("islz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%islz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%islz(iop)%evalz(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine islz

  subroutine isrx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opx ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("isrx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%isrx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%isrx(iop)%evalx(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine isrx

  subroutine isry(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opy ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("isry")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%isry(2)%ar) ) iop = 2
    endif
    CALL compact_ops%isry(iop)%evaly(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine isry

  subroutine isrz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opz ) then
      iv = v  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("isrz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%isrz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%isrz(iop)%evalz(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine isrz

  subroutine icfx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opx ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("isfx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imcfx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imcfx(iop)%evalcfx(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine icfx

  subroutine icfy(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opy ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("icfy")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imcfy(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imcfy(iop)%evalcfy(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine icfy

  subroutine icfz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opz ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("icfz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imcfz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imcfz(iop)%evalcfz(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine icfz

  subroutine ifcx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! boundary values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opx ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("ifcx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imfcx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imfcx(iop)%evalfcx(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine ifcx

  subroutine ifcy(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! boundary values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opy ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("ifcy")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imfcy(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imfcy(iop)%evalfcy(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine ifcy

  subroutine ifcz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! boundary values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    ! perform operation?
    if( compact_ops%control%null_opz ) then
!      iv = zero  ! ?
      return
    endif
    ! symmetry option
    call nvtxStartRange("ifcz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%imfcz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%imfcz(iop)%evalfcz(iv,v,vb1,vb2,iv1,iv2)
    call nvtxEndRange()
  end subroutine ifcz

  subroutine amrcfx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opx ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrcfx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrcfx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrcfx(iop)%evalx(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrcfx

  subroutine amrfcx(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opx ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrfcx")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrfcx(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrfcx(iop)%evalx(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrfcx

  subroutine amrcfy(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opy ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrcfy")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrcfy(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrcfy(iop)%evaly(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrcfy

  subroutine amrfcy(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opy ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrfcy")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrfcy(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrfcy(iop)%evaly(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrfcy

  subroutine amrcfz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opz ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrcfz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrcfz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrcfz(iop)%evalz(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrcfz

  subroutine amrfcz(v,iv,bc,vb1,vb2,iv1,iv2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: iv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: iv1,iv2 ! boundary values
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k
    if( compact_ops%control%null_opz ) then
      iv = v
      return
    endif
    ! symmetry option
    call nvtxStartRange("amrfcz")
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. allocated(compact_ops%amrfcz(2)%ar) ) iop = 2
    endif
    CALL compact_ops%amrfcz(iop)%evalz(iv,v,vb1,vb2)
    call nvtxEndRange()
  end subroutine amrfcz

end module LES_compact_operators

