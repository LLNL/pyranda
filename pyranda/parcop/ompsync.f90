module LES_ompsync
  USE iso_c_binding
  
  implicit none
  integer :: sync_var

  type fexl_pool_var
     REAL(c_double), dimension(:,:,:), allocatable :: omp_pool
     INTEGER(c_int) :: active = 0
     INTEGER(c_int) :: allocated = 0
     INTEGER(c_int) :: ax,ay,az

   contains
     procedure :: allot => allocate_pool_var
  end type fexl_pool_var


  integer(c_int) :: maxPool = 30
  type(fexl_pool_var), target :: fexlPool(1:30)
  integer(c_int) :: poolCount = 1
  
  ! subroutines
  contains
  
  subroutine allocate_pool_var(FPvar,nx,ny,nz,pc)
    implicit none
    class(fexl_pool_var), intent(out) :: FPvar
    integer(c_int),       intent(in)  :: nx,ny,nz,pc

    allocate( FPvar%omp_pool( nx,ny,nz) )
    FPvar%ax = nx
    FPvar%ay = ny
    FPvar%az = nz
    FPvar%active = 1
    FPvar%allocated = 1


    !! MAP DATA TO DEVICE HERE

    !$omp target enter data map(alloc: fexlPool(pc)%omp_pool )

    !! MAP DATA TO DEVICE HERE


    !!! $omp target enter data map(alloc: FPvar%omp_pool ) ! THis version doesnt work (btw)
    
  end subroutine allocate_pool_var

  subroutine free_pool_var(FPvar,nx,ny,nz)
    implicit none
    class(fexl_pool_var), intent(out) :: FPvar
    integer(c_int),       intent(in)  :: nx,ny,nz

    deallocate( FPvar%omp_pool )
    FPvar%active = 0
    
  end subroutine free_pool_var

  subroutine fexlPool_setPoolDepth(depth,nx,ny,nz)
    implicit none
    ! For given depth, grow the depth
    integer(c_int), intent(in) :: depth,nx,ny,nz
    integer(c_int) :: dd

    do dd=1,depth

       if ( fexlPool(dd)%allocated == 0 ) then
          
          CALL fexlPool(dd)%allot(nx,ny,nz,dd)
          fexlPool(dd)%active = 0
          
       end if
    end do
        
  end subroutine fexlPool_setPoolDepth
  
  function fexlPool_get(nx,ny,nz)
    implicit none
    real(c_double), dimension(:,:,:), pointer :: fexlPool_get
    integer(c_int), intent(in) :: nx,ny,nz
    integer(c_int) :: next

    ! Returns the next avaiable pointer
    next = poolCount + 1

    ! Need to add logic checking
    !if ( fexlPool(nextCount)%active == 1 ) then

    if ( fexlPool(next)%allocated == 0 ) then
       CALL fexlPool(next)%allot(nx,ny,nz,next)
    end if
    
    fexlPool_get => fexlPool(next)%omp_pool
    fexlPool%active = 1
    poolCount = next

    !print*,poolCount
    
  end function fexlPool_get

 
  function fexlPool_free(nitems)
    implicit none
    integer(c_int) :: fexlPool_free
    integer(c_int), intent(in)  :: nitems
    integer :: i,old

    do i = 1 , nitems

       fexlPool(poolCount)%active = 0       
       poolCount = poolCount - 1
       
    end do

    fexlPool_free = 1
       
  end function fexlPool_free
  
  
end module LES_ompsync
