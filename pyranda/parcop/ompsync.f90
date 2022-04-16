module LES_ompsync

  USE LES_objects, ONLY : pool_ptr,pool_data
  USE LES_pool
  implicit none
  integer :: sync_var

  

contains


         
! omp-pool maps
         
  function fexlPool_get() result(varptr)
    implicit none
    real(c_double), dimension(:,:,:  ), pointer :: varptr
    integer(c_int) :: next
 
    ! Returns the next avaiable pointer
    next = pool_ptr%fexlPoolCount + 1
 
    if ( next >= pool_ptr%maxPool ) then
       print*,"Max pool size reached"
       STOP
    end if
 
    print*,"fexlPool === activating ",next
     
    if ( pool_ptr%fexlPool(next)%allocated == 0 ) then
 
       !!------ ALLOCATE AND MAP --------
       allocate( pool_ptr%fexlPool(next)%omp_pool( pool_ptr%nx,pool_ptr%ny,pool_ptr%nz  ) )
       pool_ptr%fexlPool(next)%active = 1
       pool_ptr%fexlPool(next)%allocated = 1
       
       !! MAP DATA TO DEVICE HERE
       !$omp target enter data map(alloc: pool_ptr%fexlPool(next)%omp_pool )
       !! MAP DATA TO DEVICE HERE
       !!--------------------------------------
        
    print*,"       fexlPool === allocating and mapping ",next
 
 
    end if
     
    varptr => pool_ptr%fexlPool(next)%omp_pool
    pool_ptr%fexlPool(next)%active = 1
    pool_ptr%fexlPoolCount = next
     
  end function fexlPool_get
 
  function fexlPool_free(nitems,dealloc) result(error)
    implicit none
    integer(c_int) :: error
    integer(c_int), intent(in)  :: nitems             ! If (-1) will loop through all active
    integer(c_int), intent(in), optional :: dealloc   ! 1-deletes, 0-just frees
    integer :: i,old,iter,alloced,next
     
    !$omp barrier
 
    ! If nitems is -1, loop through all of the items
    if (nitems == -1) then
       iter = pool_ptr%maxPool
       pool_ptr%fexlPoolCount = pool_ptr%maxPool
    else
       iter = nitems
    end if
 
    do i = 1 , iter
 
       ! Deactivate this variable
       next = pool_ptr%fexlPoolCount
       pool_ptr%fexlPool(next)%active = 0
       alloced = pool_ptr%fexlPool(next)%allocated
        
       print*,"fexlPool === deactivate ",next
 
       ! If requested (and allocated), deallocate
       if ( present(dealloc) .and. dealloc==1 .and. alloced==1 ) then
 
          !!------ DEALLOCATE AND UNMAP --------
          !! MAP DATA TO DEVICE HERE
          !$omp target exit data map(release: pool_ptr%fexlPool(next)%omp_pool )
          !! MAP DATA TO DEVICE HERE
 
          if (allocated(pool_ptr%fexlPool(next)%omp_pool)) then
             deallocate(pool_ptr%fexlPool(next)%omp_pool )
          end if
           
          pool_ptr%fexlPool(next)%allocated = 0       
          pool_ptr%fexlPool(next)%active = 0
          !!--------------------------------------
           
          print*,"         fexlPool === deallocate & unmap ",next
 
       end if
 
       ! Increment counter
       pool_ptr%fexlPoolCount = pool_ptr%fexlPoolCount - 1
        
    end do
 
    error = 0
 
  end function fexlPool_free

         
  function fexlPoolRK4_get() result(varptr)
    implicit none
    real(c_double), dimension(:,:,: ,: ), pointer :: varptr
    integer(c_int) :: next
 
    ! Returns the next avaiable pointer
    next = pool_ptr%fexlPoolRK4Count + 1
 
    if ( next >= pool_ptr%maxPool ) then
       print*,"Max pool size reached"
       STOP
    end if
 
    print*,"fexlPoolRK4 === activating ",next
     
    if ( pool_ptr%fexlPoolRK4(next)%allocated == 0 ) then
 
       !!------ ALLOCATE AND MAP --------
       allocate( pool_ptr%fexlPoolRK4(next)%omp_pool( pool_ptr%nx,pool_ptr%ny,pool_ptr%nz ,pool_ptr%nRK4 ) )
       pool_ptr%fexlPoolRK4(next)%active = 1
       pool_ptr%fexlPoolRK4(next)%allocated = 1
       
       !! MAP DATA TO DEVICE HERE
       !$omp target enter data map(alloc: pool_ptr%fexlPoolRK4(next)%omp_pool )
       !! MAP DATA TO DEVICE HERE
       !!--------------------------------------
        
    print*,"       fexlPoolRK4 === allocating and mapping ",next
 
 
    end if
     
    varptr => pool_ptr%fexlPoolRK4(next)%omp_pool
    pool_ptr%fexlPoolRK4(next)%active = 1
    pool_ptr%fexlPoolRK4Count = next
     
  end function fexlPoolRK4_get
 
  function fexlPoolRK4_free(nitems,dealloc) result(error)
    implicit none
    integer(c_int) :: error
    integer(c_int), intent(in)  :: nitems             ! If (-1) will loop through all active
    integer(c_int), intent(in), optional :: dealloc   ! 1-deletes, 0-just frees
    integer :: i,old,iter,alloced,next
     
    !$omp barrier
 
    ! If nitems is -1, loop through all of the items
    if (nitems == -1) then
       iter = pool_ptr%maxPool
       pool_ptr%fexlPoolRK4Count = pool_ptr%maxPool
    else
       iter = nitems
    end if
 
    do i = 1 , iter
 
       ! Deactivate this variable
       next = pool_ptr%fexlPoolRK4Count
       pool_ptr%fexlPoolRK4(next)%active = 0
       alloced = pool_ptr%fexlPoolRK4(next)%allocated
        
       print*,"fexlPoolRK4 === deactivate ",next
 
       ! If requested (and allocated), deallocate
       if ( present(dealloc) .and. dealloc==1 .and. alloced==1 ) then
 
          !!------ DEALLOCATE AND UNMAP --------
          !! MAP DATA TO DEVICE HERE
          !$omp target exit data map(release: pool_ptr%fexlPoolRK4(next)%omp_pool )
          !! MAP DATA TO DEVICE HERE
 
          if (allocated(pool_ptr%fexlPoolRK4(next)%omp_pool)) then
             deallocate(pool_ptr%fexlPoolRK4(next)%omp_pool )
          end if
           
          pool_ptr%fexlPoolRK4(next)%allocated = 0       
          pool_ptr%fexlPoolRK4(next)%active = 0
          !!--------------------------------------
           
          print*,"         fexlPoolRK4 === deallocate & unmap ",next
 
       end if
 
       ! Increment counter
       pool_ptr%fexlPoolRK4Count = pool_ptr%fexlPoolRK4Count - 1
        
    end do
 
    error = 0
 
  end function fexlPoolRK4_free

         
  function fexlPoolSpecies_get() result(varptr)
    implicit none
    real(c_double), dimension(:,:,: ,: ), pointer :: varptr
    integer(c_int) :: next
 
    ! Returns the next avaiable pointer
    next = pool_ptr%fexlPoolSpeciesCount + 1
 
    if ( next >= pool_ptr%maxPool ) then
       print*,"Max pool size reached"
       STOP
    end if
 
    print*,"fexlPoolSpecies === activating ",next
     
    if ( pool_ptr%fexlPoolSpecies(next)%allocated == 0 ) then
 
       !!------ ALLOCATE AND MAP --------
       allocate( pool_ptr%fexlPoolSpecies(next)%omp_pool( pool_ptr%nx,pool_ptr%ny,pool_ptr%nz ,pool_ptr%nSpecies ) )
       pool_ptr%fexlPoolSpecies(next)%active = 1
       pool_ptr%fexlPoolSpecies(next)%allocated = 1
       
       !! MAP DATA TO DEVICE HERE
       !$omp target enter data map(alloc: pool_ptr%fexlPoolSpecies(next)%omp_pool )
       !! MAP DATA TO DEVICE HERE
       !!--------------------------------------
        
    print*,"       fexlPoolSpecies === allocating and mapping ",next
 
 
    end if
     
    varptr => pool_ptr%fexlPoolSpecies(next)%omp_pool
    pool_ptr%fexlPoolSpecies(next)%active = 1
    pool_ptr%fexlPoolSpeciesCount = next
     
  end function fexlPoolSpecies_get
 
  function fexlPoolSpecies_free(nitems,dealloc) result(error)
    implicit none
    integer(c_int) :: error
    integer(c_int), intent(in)  :: nitems             ! If (-1) will loop through all active
    integer(c_int), intent(in), optional :: dealloc   ! 1-deletes, 0-just frees
    integer :: i,old,iter,alloced,next
     
    !$omp barrier
 
    ! If nitems is -1, loop through all of the items
    if (nitems == -1) then
       iter = pool_ptr%maxPool
       pool_ptr%fexlPoolSpeciesCount = pool_ptr%maxPool
    else
       iter = nitems
    end if
 
    do i = 1 , iter
 
       ! Deactivate this variable
       next = pool_ptr%fexlPoolSpeciesCount
       pool_ptr%fexlPoolSpecies(next)%active = 0
       alloced = pool_ptr%fexlPoolSpecies(next)%allocated
        
       print*,"fexlPoolSpecies === deactivate ",next
 
       ! If requested (and allocated), deallocate
       if ( present(dealloc) .and. dealloc==1 .and. alloced==1 ) then
 
          !!------ DEALLOCATE AND UNMAP --------
          !! MAP DATA TO DEVICE HERE
          !$omp target exit data map(release: pool_ptr%fexlPoolSpecies(next)%omp_pool )
          !! MAP DATA TO DEVICE HERE
 
          if (allocated(pool_ptr%fexlPoolSpecies(next)%omp_pool)) then
             deallocate(pool_ptr%fexlPoolSpecies(next)%omp_pool )
          end if
           
          pool_ptr%fexlPoolSpecies(next)%allocated = 0       
          pool_ptr%fexlPoolSpecies(next)%active = 0
          !!--------------------------------------
           
          print*,"         fexlPoolSpecies === deallocate & unmap ",next
 
       end if
 
       ! Increment counter
       pool_ptr%fexlPoolSpeciesCount = pool_ptr%fexlPoolSpeciesCount - 1
        
    end do
 
    error = 0
 
  end function fexlPoolSpecies_free

         
  function fexlPoolVector_get() result(varptr)
    implicit none
    real(c_double), dimension(:,:,: ,: ), pointer :: varptr
    integer(c_int) :: next
 
    ! Returns the next avaiable pointer
    next = pool_ptr%fexlPoolVectorCount + 1
 
    if ( next >= pool_ptr%maxPool ) then
       print*,"Max pool size reached"
       STOP
    end if
 
    print*,"fexlPoolVector === activating ",next
     
    if ( pool_ptr%fexlPoolVector(next)%allocated == 0 ) then
 
       !!------ ALLOCATE AND MAP --------
       allocate( pool_ptr%fexlPoolVector(next)%omp_pool( pool_ptr%nx,pool_ptr%ny,pool_ptr%nz ,3 ) )
       pool_ptr%fexlPoolVector(next)%active = 1
       pool_ptr%fexlPoolVector(next)%allocated = 1
       
       !! MAP DATA TO DEVICE HERE
       !$omp target enter data map(alloc: pool_ptr%fexlPoolVector(next)%omp_pool )
       !! MAP DATA TO DEVICE HERE
       !!--------------------------------------
        
    print*,"       fexlPoolVector === allocating and mapping ",next
 
 
    end if
     
    varptr => pool_ptr%fexlPoolVector(next)%omp_pool
    pool_ptr%fexlPoolVector(next)%active = 1
    pool_ptr%fexlPoolVectorCount = next
     
  end function fexlPoolVector_get
 
  function fexlPoolVector_free(nitems,dealloc) result(error)
    implicit none
    integer(c_int) :: error
    integer(c_int), intent(in)  :: nitems             ! If (-1) will loop through all active
    integer(c_int), intent(in), optional :: dealloc   ! 1-deletes, 0-just frees
    integer :: i,old,iter,alloced,next
     
    !$omp barrier
 
    ! If nitems is -1, loop through all of the items
    if (nitems == -1) then
       iter = pool_ptr%maxPool
       pool_ptr%fexlPoolVectorCount = pool_ptr%maxPool
    else
       iter = nitems
    end if
 
    do i = 1 , iter
 
       ! Deactivate this variable
       next = pool_ptr%fexlPoolVectorCount
       pool_ptr%fexlPoolVector(next)%active = 0
       alloced = pool_ptr%fexlPoolVector(next)%allocated
        
       print*,"fexlPoolVector === deactivate ",next
 
       ! If requested (and allocated), deallocate
       if ( present(dealloc) .and. dealloc==1 .and. alloced==1 ) then
 
          !!------ DEALLOCATE AND UNMAP --------
          !! MAP DATA TO DEVICE HERE
          !$omp target exit data map(release: pool_ptr%fexlPoolVector(next)%omp_pool )
          !! MAP DATA TO DEVICE HERE
 
          if (allocated(pool_ptr%fexlPoolVector(next)%omp_pool)) then
             deallocate(pool_ptr%fexlPoolVector(next)%omp_pool )
          end if
           
          pool_ptr%fexlPoolVector(next)%allocated = 0       
          pool_ptr%fexlPoolVector(next)%active = 0
          !!--------------------------------------
           
          print*,"         fexlPoolVector === deallocate & unmap ",next
 
       end if
 
       ! Increment counter
       pool_ptr%fexlPoolVectorCount = pool_ptr%fexlPoolVectorCount - 1
        
    end do
 
    error = 0
 
  end function fexlPoolVector_free




        SUBROUTINE ompdata_set_pool_rk4(nrk4)
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: nrk4
          
          pool_ptr%nRK4 = nrk4
          
        END SUBROUTINE ompdata_set_pool_rk4


        
        subroutine ompdata_removePools(patch,level)
          integer, intent(in) :: patch,level
          integer :: error


          pool_ptr => pool_data(patch,level)
          
          error = fexlPool_free(-1,1)
          error = fexlPoolVector_free(-1,1)
          error = fexlPoolSpecies_free(-1,1)
          error = fexlPoolRK4_free(-1,1)
          
        end subroutine ompdata_removePools

        



end module LES_ompsync
