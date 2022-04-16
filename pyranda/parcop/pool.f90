!===================================================================================================
MODULE LES_pool
!===================================================================================================
  USE iso_c_binding
  USE LES_patch, ONLY : patch_type
  
  implicit none
  !integer :: sync_var
  integer(c_int) :: pool_error
  
  type pool_var
     REAL(c_double), dimension(:,:,:), allocatable :: omp_pool
     INTEGER(c_int) :: active = 0
     INTEGER(c_int) :: allocated = 0
  end type pool_var

  type pool_var4d
     REAL(c_double), dimension(:,:,:,:), allocatable :: omp_pool
     INTEGER(c_int) :: active = 0
     INTEGER(c_int) :: allocated = 0
  end type pool_var4d

    
  type pool_type

     ! Global sizes for set pools
     integer(c_int)   :: maxPool = 50
     type(pool_var)   :: fexlPool(1:50)
     type(pool_var4d) :: fexlPoolVector(1:50)
     type(pool_var4d) :: fexlPoolSpecies(1:50)
     type(pool_var4d) :: fexlPoolRK4(1:50)
     integer(c_int) :: fexlPoolCount         = 0
     integer(c_int) :: fexlPoolVectorCount   = 0
     integer(c_int) :: fexlPoolSpeciesCount  = 0
     integer(c_int) :: fexlPoolRK4Count      = 0
     integer(c_int) :: nx, ny, nz
     integer(c_int) :: nSpecies,nRK4


     contains
       procedure :: setup => setup_pool
       procedure :: remove => remove_pool     
  end type pool_type
!===================================================================================================
  contains
!===================================================================================================
    ! PASS
    !   access functions are in ompdata.py.f90
    subroutine setup_pool(pool_data,patch_data,ns)
      class(pool_type), intent(inout) :: pool_data
      class(patch_type), intent(in)   :: patch_data
      integer, intent(in)             :: ns

      
      CALL remove_pool(pool_data)
      
      ! Copy over array sizes
      pool_data%nx = patch_data%ax
      pool_data%ny = patch_data%ay
      pool_data%nz = patch_data%az
      pool_data%nspecies = ns
      
    end subroutine setup_pool


    subroutine remove_pool(pool_data)
      class(pool_type), intent(inout) :: pool_data
      integer :: ip
      
      do ip=1,pool_data%maxPool
         
         if (allocated(pool_data%fexlPool(ip)%omp_pool)) then
            deallocate(pool_data%fexlPool(ip)%omp_pool)
         end if
         if (allocated(pool_data%fexlPoolVector(ip)%omp_pool)) then
            deallocate(pool_data%fexlPoolVector(ip)%omp_pool)
         end if
         if (allocated(pool_data%fexlPoolSpecies(ip)%omp_pool)) then
            deallocate(pool_data%fexlPoolSpecies(ip)%omp_pool)
         end if
         if (allocated(pool_data%fexlPoolRK4(ip)%omp_pool)) then
            deallocate(pool_data%fexlPoolRK4(ip)%omp_pool)
         end if

         pool_data%fexlPool(ip)%active = 0
         pool_data%fexlPool(ip)%allocated = 0

         pool_data%fexlPoolVector(ip)%active = 0
         pool_data%fexlPoolVector(ip)%allocated = 0

         pool_data%fexlPoolSpecies(ip)%active = 0
         pool_data%fexlPoolSpecies(ip)%allocated = 0

         pool_data%fexlPoolRK4(ip)%active = 0
         pool_data%fexlPoolRK4(ip)%allocated = 0
         
      end do

      pool_data%fexlPoolCount = 0
      pool_data%fexlPoolVectorCount = 0
      pool_data%fexlPoolSpeciesCount = 0
      pool_data%fexlPoolRK4Count = 0
      
    end subroutine remove_pool

    
end module LES_pool
