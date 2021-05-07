!===================================================================================================
 MODULE LES_pentadiagonal ! pentadiagonal matrix solvers
  USE iso_c_binding
  !use LES_input, only : gpu_kernel
  USE LES_ompsync
  interface btrid_block4_lus
    module procedure btrid_block4_lus_al
  end interface

  interface ptrid_block4_lus
    module procedure ptrid_block4_lus_al
  end interface

  integer :: gpu_kernel = 1
  
  contains

  subroutine bpentLUD1(c,n)
    implicit none
    integer(c_int), intent(in) :: n
    real(kind=c_double), dimension(n,5), intent(inout) :: c
    integer(c_int) :: i
    ! Forward elimination
    do i=1,n-2
      ! Eliminate b(i+1)
      c(i+1,2) = c(i+1,2)/c(i,3)
      c(i+1,3) = c(i+1,3) - c(i,4)*c(i+1,2)
      c(i+1,4) = c(i+1,4) - c(i,5)*c(i+1,2)
      ! Eliminate a(i+2)
      c(i+2,1) = c(i+2,1)/c(i,3)
      c(i+2,2) = c(i+2,2) - c(i,4)*c(i+2,1)
      c(i+2,3) = c(i+2,3) - c(i,5)*c(i+2,1)
    end do
    ! Eliminate b(n)
    c(n,2) = c(n,2)/c(n-1,3)
    c(n,3) = c(n,3) - c(n-1,4)*c(n,2)
    ! Pre-divide diagonal
    c(:,3) = 1.0_c_double/c(:,3)
  end subroutine bpentLUD1


  subroutine bpentLUS1(c, r, n)
    implicit none
    integer(c_int), intent(in) :: n
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n), intent(inout) :: r
    integer(c_int) :: i
    ! Forward substitution
    do i=1,n-2
      r(i+1) = r(i+1) - r(i)*c(i+1,2)
      r(i+2) = r(i+2) - r(i)*c(i+2,1)
    end do
    r(n) = r(n) - r(n-1)*c(n,2)
    ! Backward substitution
    ! Diagonal has been pre-divided
    r(n) = r(n)*c(n,3)
    r(n-1) = (r(n-1) - c(n-1,4)*r(n))*c(n-1,3)
    do i=n-2,1,-1
      r(i) = (r(i) - c(i,4)*r(i+1) - c(i,5)*r(i+2))*c(i,3)
    end do
  end subroutine bpentLUS1


  SUBROUTINE ppentLUD1(c,N)
      IMPLICIT NONE
      integer(c_int), INTENT(IN) :: N
      real(kind=c_double), dimension(N,9), INTENT(INOUT) :: c
      real(kind=c_double), parameter :: one=1.0_c_double
      integer(c_int) :: K
  
      c(1,8)=c(1,1)
      c(1,9)=c(1,2)
      c(1,3)=one/c(1,3)
      c(1,6)=c(N-1,5)*c(1,3)
      c(1,7)=c(N,4)*c(1,3)

      c(2,2)=c(2,2)*c(1,3)
      c(2,3)=c(2,3)-c(2,2)*c(1,4)
      c(2,4)=c(2,4)-c(2,2)*c(1,5)
      c(2,8)=-c(2,2)*c(1,8)
      c(2,9)=c(2,1)-c(2,2)*c(1,9)
      c(2,3)=one/c(2,3)
      c(2,6)=-c(1,6)*c(1,4)*c(2,3)
      c(2,7)=(c(N,5)-c(1,7)*c(1,4))*c(2,3)

      DO K=3,N-4
        c(K,1)=c(K,1)*c(K-2,3)
        c(K,2)=(c(K,2)-c(K,1)*c(K-2,4))*c(K-1,3)
        c(K,3)=c(K,3)-(c(K,2)*c(K-1,4)+c(K,1)*c(K-2,5))
        c(K,4)=c(K,4)-c(K,2)*c(K-1,5)
        c(K,8)=-(c(K,2)*c(K-1,8)+c(K,1)*c(K-2,8))
        c(K,9)=-(c(K,2)*c(K-1,9)+c(K,1)*c(K-2,9))
        c(K,3)=one/c(K,3)
        c(K,6)=-(c(K-1,6)*c(K-1,4)+c(K-2,6)*c(K-2,5))*c(K,3)
        c(K,7)=-(c(K-1,7)*c(K-1,4)+c(K-2,7)*c(K-2,5))*c(K,3)
      END DO

      c(N-3,1)=c(N-3,1)*c(N-5,3)
      c(N-3,2)=(c(N-3,2)-c(N-3,1)*c(N-5,4))*c(N-4,3)
      c(N-3,3)=c(N-3,3)-(c(N-3,2)*c(N-4,4)+c(N-3,1)*c(N-5,5))
      c(N-3,4)=c(N-3,4)-c(N-3,2)*c(N-4,5)
      c(N-3,8)=c(N-3,5)-(c(N-3,2)*c(N-4,8)+c(N-3,1)*c(N-5,8))
      c(N-3,9)=-(c(N-3,2)*c(N-4,9)+c(N-3,1)*c(N-5,9))
      c(N-3,3)=one/c(N-3,3)
      c(N-3,6)=(c(N-1,1)-(c(N-4,6)*c(N-4,4)+c(N-5,6)*c(N-5,5)))*c(N-3,3)
      c(N-3,7)=-(c(N-4,7)*c(N-4,4)+c(N-5,7)*c(N-5,5))*c(N-3,3)

      c(N-2,1)=c(N-2,1)*c(N-4,3)
      c(N-2,2)=(c(N-2,2)-c(N-2,1)*c(N-4,4))*c(N-3,3)
      c(N-2,3)=c(N-2,3)-(c(N-2,2)*c(N-3,4)+c(N-2,1)*c(N-4,5))
      c(N-2,8)=c(N-2,4)-(c(N-2,2)*c(N-3,8)+c(N-2,1)*c(N-4,8))
      c(N-2,9)=c(N-2,5)-(c(N-2,2)*c(N-3,9)+c(N-2,1)*c(N-4,9))
      c(N-2,3)=one/c(N-2,3)
      c(N-2,6)=(c(N-1,2)-c(N-3,6)*c(N-3,4)-c(N-4,6)*c(N-4,5))*c(N-2,3)
      c(N-2,7)=(c(N,1)-c(N-3,7)*c(N-3,4)-c(N-4,7)*c(N-4,5))*c(N-2,3)

      DO K=1,N-2
        c(N-1,3)=c(N-1,3)-c(K,6)*c(K,8)
        c(N-1,4)=c(N-1,4)-c(K,6)*c(K,9)
      END DO
      c(N-1,9)=c(N-1,4)
      DO K=1,N-2
        c(N,2)=c(N,2)-c(K,7)*c(K,8)
      END DO
      c(N-1,3)=one/c(N-1,3)
      c(N-1,7)=c(N,2)*c(N-1,3)
      DO K=1,N-1
        c(N,3)=c(N,3)-c(K,7)*c(K,9)
      END DO
      c(N,3)=one/c(N,3)
  END SUBROUTINE ppentLUD1


  SUBROUTINE ppentLUS1(c,r,n)
      IMPLICIT NONE
      integer(c_int), intent(in) :: n
      real(kind=c_double), DIMENSION(n,9), INTENT(IN) :: c
      real(kind=c_double), DIMENSION(n), INTENT(INOUT) :: r
      integer(c_int) :: i
        r(2)=r(2)-c(2,2)*r(1)
        DO i=3,N-2
          r(i)=r(i)-(c(i,2)*r(i-1)+c(i,1)*r(i-2))
        END DO
        DO i=1,N-2
          r(N-1)=r(N-1)-c(i,6)*r(i)
        END DO
        DO i=1,N-1
          r(N)=r(N)-c(i,7)*r(i)
        END DO
        r(N)=r(N)*c(N,3)
        r(N-1)=(r(N-1)-c(N-1,9)*r(N))*c(N-1,3)
        r(N-2)=(r(N-2)-c(N-2,8)*r(N-1)-c(N-2,9)*r(N))*c(N-2,3)
        r(N-3)=(r(N-3)-(c(N-3,4)*r(N-2)+c(N-3,8)*r(N-1)+c(N-3,9)*r(N)))*c(N-3,3)
        DO i=N-4,1,-1
          r(i)=(r(i)-(c(i,4)*r(i+1)+c(i,5)*r(i+2)+c(i,8)*r(n-1)+c(i,9)*r(n)))*c(i,3)
        END DO
  END SUBROUTINE ppentLUS1


  subroutine btrid_block4_lud( abc, n )
    USE LES_blockmath
    implicit none
    integer(c_int), intent(in) :: n
    real(kind=c_double), dimension(4,4,3,n), intent(inout) :: abc
    type(block4), dimension(n) :: a,b,c
    integer(c_int) :: j
    do j=1,n
      a(j) = abc(:,:,1,j)
      b(j) = abc(:,:,2,j)
      c(j) = abc(:,:,3,j)
    end do
    b(1) = .inv.b(1)
    abc(:,:,2,1) = b(1)%array
    do j=2,n
      a(j) = a(j)*b(j-1)
      b(j) = .inv.( b(j) - a(j)*c(j-1) )
      abc(:,:,1,j) = a(j)%array
      abc(:,:,2,j) = b(j)%array
    end do
  end subroutine btrid_block4_lud


! note: a(:,3:4,3,:) = 0
  subroutine btrid_block4_lus_al( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,3,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4) :: tmp
    integer(c_int) :: i,j,k
    !$omp target teams if(gpu_kernel==1) !$fexl-async.sync_var.
    !$omp distribute parallel do collapse(2)
    do j=1,n2
     do i=1,n1
      do k=2,n
        r(:,i,j,k) = r(:,i,j,k) - a(:,1,1,k)*r(1,i,j,k-1) - a(:,2,1,k)*r(2,i,j,k-1) &
                                - a(:,3,1,k)*r(3,i,j,k-1) - a(:,4,1,k)*r(4,i,j,k-1)
      end do
     end do
    end do
    !$omp end distribute parallel do
    !$omp distribute parallel do collapse(2) private(tmp)
    do j=1,n2
     do i=1,n1
	  	 tmp(:) = r(:,i,j,n)
       r(:,i,j,n) = a(:,1,2,n)*tmp(1) + a(:,2,2,n)*tmp(2) &
       	          + a(:,3,2,n)*tmp(3) + a(:,4,2,n)*tmp(4)
     end do
    end do
    !$omp end distribute parallel do
    !$omp distribute parallel do collapse(2) private(tmp)
     do j=1,n2
      do i=1,n1
      do k=n-1,1,-1
	   		tmp(:) = r(:,i,j,k) - a(:,1,3,k)*r(1,i,j,k+1) - a(:,2,3,k)*r(2,i,j,k+1)
       	r(:,i,j,k) = a(:,1,2,k)*tmp(1) + a(:,2,2,k)*tmp(2) &
                   + a(:,3,2,k)*tmp(3) + a(:,4,2,k)*tmp(4)
      end do
     end do
    end do
    !$omp end distribute parallel do
    !$omp end target teams
  end subroutine btrid_block4_lus_al

! note: a(:,3:4,3,:) = 0
  subroutine btrid_block4_lus_as( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,3,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4) :: tmp
    integer(c_int) :: i,j,k,l
    do k=2,n
      do j=1,n2
      do i=1,n1
        do l=1,4
          r(l,i,j,k) = r(l,i,j,k) - sum(a(l,:,1,k)*r(:,i,j,k-1))
        end do
      end do
      end do
    end do
    do j=1,n2
    do i=1,n1
      do l=1,4
        tmp(l) = sum(a(l,:,2,n)*r(:,i,j,n))
      end do
      r(:,i,j,n) = tmp
    end do
    end do
    do k=n-1,1,-1
      do j=1,n2
      do i=1,n1
      do l=1,4
        tmp(l) = sum( a(l,:,2,k)*( r(:,i,j,k) - a(:,1,3,k)*r(1,i,j,k+1) - a(:,2,3,k)*r(2,i,j,k+1) ) )
      end do
      r(:,i,j,k) = tmp
      end do
      end do
    end do
  end subroutine btrid_block4_lus_as

! note: a(:,3:4,3,:) = 0
  subroutine btrid_block4_lus_mv( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,3,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    integer(c_int) :: i,j,k
    do k=2,n
      do j=1,n2
      do i=1,n1
        r(:,i,j,k) = r(:,i,j,k) - matmul(a(:,:,1,k),r(:,i,j,k-1))
      end do
      end do
    end do
    do j=1,n2
    do i=1,n1
      r(:,i,j,n) = matmul(a(:,:,2,n),r(:,i,j,n))
    end do
    end do
    do k=n-1,1,-1
      do j=1,n2
      do i=1,n1
        r(:,i,j,k) = matmul( a(:,:,2,k), r(:,i,j,k) - matmul(a(:,1:2,3,k),r(1:2,i,j,k+1)) )
      end do
      end do
    end do
  end subroutine btrid_block4_lus_mv


! use dgemm-able version
! note: a(:,3:4,3,:) = 0
  subroutine btrid_block4_lus_dgemm( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,3,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4,n1) :: tmp1,tmp2
    real(kind=c_double), dimension(4,4) :: c
    integer(c_int) :: i,j,k
	c = 0.0_8
    do k=2,n
      do j=1,n2
        tmp1 = r(:,:,j,k-1)
        tmp2 = matmul(a(:,:,1,k),tmp1)
        r(:,:,j,k) = r(:,:,j,k) - tmp2
      end do
    end do
    do j=1,n2
      tmp1 = r(:,:,j,n)
      tmp2 = matmul(a(:,:,2,n),tmp1)
      r(:,:,j,n) = tmp2
    end do
    do k=n-1,1,-1
	  c(:,1:2) = a(:,1:2,3,k)
      do j=1,n2
        tmp1 = r(:,:,j,k+1)
        tmp2 = matmul(c,tmp1)
        tmp1 = r(:,:,j,k) - tmp2
        tmp2 = matmul(a(:,:,2,k),tmp1)
        r(:,:,j,k) = tmp2
      end do
    end do
  end subroutine btrid_block4_lus_dgemm


! note: a(:,3:4,3,:) = 0
  subroutine btrid_block4_lus_mm( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,3,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    integer(c_int) :: i,j,k
    do k=2,n
      do j=1,n2
        r(:,:,j,k) = r(:,:,j,k) - matmul(a(:,:,1,k),r(:,:,j,k-1))
      end do
    end do
    do j=1,n2
      r(:,:,j,n) = matmul(a(:,:,2,n),r(:,:,j,n))
    end do
    do k=n-1,1,-1
      do j=1,n2
        r(:,:,j,k) = matmul( a(:,:,2,k), r(:,:,j,k) - matmul(a(:,1:2,3,k),r(1:2,:,j,k+1)) )
      end do
    end do
  end subroutine btrid_block4_lus_mm


  subroutine ptrid_block4_lud( abc, n )
    USE LES_blockmath
    implicit none
    integer(c_int), intent(in) :: n
    real(kind=c_double), dimension(4,4,4,n), intent(inout) :: abc
    type(block4), dimension(n) :: a,b,c,aa,cc
    integer(c_int) :: j
    do j=1,n
      a(j) = abc(:,:,1,j)
      b(j) = abc(:,:,2,j)
      c(j) = abc(:,:,3,j)
      aa(j) = 0.0_8
      cc(j) = 0.0_8
    end do
    j = 1
    b(j) = .inv.b(j)
    aa(j) = a(j)
    cc(j) = c(n)*b(j)
    do j=2,n-1
      b(n) = b(n)-cc(j-1)*aa(j-1)
      c(n) = -cc(j-1)*c(j-1)
      a(j) = a(j)*b(j-1)
      aa(j) = -a(j)*aa(j-1)
      b(j) = .inv.( b(j) - a(j)*c(j-1) )
      cc(j) = c(n)*b(j)
    end do
    j = n
    a(j) = a(j)*b(j-1) + cc(j-1)
    b(j) = .inv.( b(j) - a(j)*(c(j-1)+aa(j-1)) )
    do j=1,n
      abc(:,:,1,j) = a(j)%array
      abc(:,:,2,j) = b(j)%array
      abc(:,1:2,4,j) = cc(j)%array(:,1:2)
      abc(:,3:4,4,j) = aa(j)%array(:,3:4)
!      abc(:,:,4,j) = aa(j)%array
!      abc(:,:,5,j) = cc(j)%array
    end do
  end subroutine ptrid_block4_lud


! note: a(:,3:4,3,:) = 0, a(:,3:4,5,:) = 0, and a(:,1:2,4,:) = 0
! a(:,1:2,5,:) packed in a(:,1:2,4,:)
  subroutine ptrid_block4_lus_al( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,4,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4) :: tmp
    integer(c_int) :: i,j,k
    !$omp target teams if(gpu_kernel==1) !$fexl-async.sync_var.
    if (n > 2) then
    !$omp distribute parallel do collapse(2)
    do j=1,n2
     do i=1,n1
       r(:,i,j,n) = r(:,i,j,n) - a(:,1,4,1)*r(1,i,j,1) - a(:,2,4,1)*r(2,i,j,1)
     end do
    end do
    !$omp end distribute parallel do
    endif
    !$omp distribute parallel do collapse(2)
    do j=1,n2
     do i=1,n1
      do k=2,n-2
        r(:,i,j,k) = r(:,i,j,k) - a(:,1,1,k)*r(1,i,j,k-1) - a(:,2,1,k)*r(2,i,j,k-1) &
                                - a(:,3,1,k)*r(3,i,j,k-1) - a(:,4,1,k)*r(4,i,j,k-1)
        r(:,i,j,n) = r(:,i,j,n) - a(:,1,4,k)*r(1,i,j,k)   - a(:,2,4,k)*r(2,i,j,k)
      end do
     end do
    end do
    !$omp end distribute parallel do
    if (n > 2) then
    !$omp distribute parallel do collapse(2)
    do j=1,n2
     do i=1,n1
       r(:,i,j,n-1) = r(:,i,j,n-1) - a(:,1,1,n-1)*r(1,i,j,n-2) - a(:,2,1,n-1)*r(2,i,j,n-2) &
                                   - a(:,3,1,n-1)*r(3,i,j,n-2) - a(:,4,1,n-1)*r(4,i,j,n-2)
     end do
    end do
    !$omp end distribute parallel do
    end if
    !$omp distribute parallel do collapse(2) private(tmp)
    do j=1,n2
     do i=1,n1
       tmp(:)     = r(:,i,j,n) - a(:,1,1,n)*r(1,i,j,n-1) - a(:,2,1,n)*r(2,i,j,n-1) &
                                 - a(:,3,1,n)*r(3,i,j,n-1) - a(:,4,1,n)*r(4,i,j,n-1)
       r(:,i,j,n) = a(:,1,2,n)*tmp(1) + a(:,2,2,n)*tmp(2) &
                  + a(:,3,2,n)*tmp(3) + a(:,4,2,n)*tmp(4)
     end do
    end do
    !$omp end distribute parallel do
    !$omp distribute parallel do collapse(2) private(tmp)
    do j=1,n2
     do i=1,n1
      do k=n-1,1,-1
	      tmp(:) = r(:,i,j,k) - a(:,1,3,k)*r(1,i,j,k+1) - a(:,2,3,k)*r(2,i,j,k+1) &
                              - a(:,3,4,k)*r(3,i,j,n)   - a(:,4,4,k)*r(4,i,j,n)
        r(:,i,j,k) = a(:,1,2,k)*tmp(1) + a(:,2,2,k)*tmp(2) &
                   + a(:,3,2,k)*tmp(3) + a(:,4,2,k)*tmp(4)
      end do
     end do
    end do
    !$omp end distribute parallel do
    !$omp end target teams
  end subroutine ptrid_block4_lus_al

  subroutine ptrid_block4_lus_as( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,4,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4) :: tmp
    integer(c_int) :: i,j,k,l
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do j=1,n2
     do i=1,n1
      do k=2,n
       do l=1,4
        r(l,i,j,k) = r(l,i,j,k) - sum(a(l,:,1,k)*r(:,i,j,k-1))
       end do ! l
      end do ! k
     end do ! i
    end do ! j
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do j=1,n2
     do i=1,n1
      do k=1,n-2
       do l=1,4
        r(l,i,j,n) = r(l,i,j,n) - a(l,1,4,k)*r(1,i,j,k) - a(l,2,4,k)*r(2,i,j,k)
       end do
      end do
     end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do collapse(2) private(tmp) if(gpu_kernel==1) !$fexl-async.sync_var.
    do j=1,n2
    do i=1,n1
      do l=1,4
        tmp(l) = sum(a(l,:,2,n)*r(:,i,j,n))
      end do
      r(:,i,j,n) = tmp
    end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do j=1,n2
     do i=1,n1
      do k=n-1,1,-1
       do l=1,4
        tmp(l) = sum(a(l,:,2,k)*( r(:,i,j,k) - a(:,1,3,k)*r(1,i,j,k+1) - a(:,2,3,k)*r(2,i,j,k+1) &
          - a(:,3,4,k)*r(3,i,j,n) - a(:,4,4,k)*r(4,i,j,n) ) )
       end do
        r(:,i,j,k) = tmp
      end do
     end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ptrid_block4_lus_as

! note: a(:,3:4,3,:) = 0, a(:,3:4,5,:) = 0, and a(:,1:2,4,:) = 0
! a(:,1:2,5,:) packed in a(:,1:2,4,:)
  subroutine ptrid_block4_lus_mv( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,4,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    integer(c_int) :: i,j,k
    do k=2,n
      do j=1,n2
      do i=1,n1
        r(:,i,j,k) = r(:,i,j,k) - matmul(a(:,:,1,k),r(:,i,j,k-1))
      end do
      end do
    end do
    do k=1,n-2
      do j=1,n2
      do i=1,n1
        r(:,i,j,n) = r(:,i,j,n) - matmul(a(:,1:2,4,k),r(1:2,i,j,k))
      end do
      end do
    end do
    do j=1,n2
    do i=1,n1
      r(:,i,j,n) = matmul(a(:,:,2,n),r(:,i,j,n))
    end do
    end do
    do k=n-1,1,-1
      do j=1,n2
      do i=1,n1
        r(:,i,j,k) = matmul( a(:,:,2,k), r(:,i,j,k)-matmul(a(:,1:2,3,k),r(1:2,i,j,k+1))-matmul(a(:,3:4,4,k),r(3:4,i,j,n)) )
      end do
      end do
    end do
  end subroutine ptrid_block4_lus_mv

! use dgemm-able version
! note: a(:,3:4,3,:) = 0, a(:,3:4,5,:) = 0, and a(:,1:2,4,:) = 0
! a(:,1:2,5,:) packed in a(:,1:2,4,:)
  subroutine ptrid_block4_lus_dgemm( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,4,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    real(kind=c_double), dimension(4,n1) :: tmp1,tmp2
    real(kind=c_double), dimension(4,4) :: c
    integer(c_int) :: i,j,k
    c = 0.0_8
    c(:,1:2) = a(:,1:2,4,1)
    do j=1,n2
      tmp1 = r(:,:,j,1)
      tmp2 = matmul(c,tmp1)
      r(:,:,j,n) = r(:,:,j,n) - tmp2
    end do
    do k=2,n-2
      c(:,1:2) = a(:,1:2,4,k)
      do j=1,n2
        tmp1 = r(:,:,j,k-1)
        tmp2 = matmul(a(:,:,1,k),tmp1)
        r(:,:,j,k) = r(:,:,j,k) - tmp2
        tmp1 = r(:,:,j,k)
        tmp2 =  matmul(c,tmp1)
        r(:,:,j,n) = r(:,:,j,n) - tmp2
      end do
    end do
    do j=1,n2
      tmp1 = r(:,:,j,n-2)
      tmp2 = matmul(a(:,:,1,n-1),tmp1)
      r(:,:,j,n-1) = r(:,:,j,n-1) - tmp2
      tmp1 = r(:,:,j,n-1)
      tmp2 = matmul(a(:,:,1,n),tmp1)
      tmp1 = r(:,:,j,n) - tmp2
      tmp2 = matmul(a(:,:,2,n),tmp1)
      r(:,:,j,n) = tmp2
    end do
    do k=n-1,1,-1
      c(:,1:2) = a(:,1:2,3,k)
      c(:,3:4) = a(:,3:4,4,k)
      do j=1,n2
        tmp1(1:2,:) = r(1:2,:,j,k+1)
        tmp1(3:4,:) = r(3:4,:,j,n)
        tmp2 = matmul(c,tmp1)
        tmp1 = r(:,:,j,k) - tmp2
        tmp2 = matmul(a(:,:,2,k),tmp1)
        r(:,:,j,k) = tmp2
      end do
    end do
  end subroutine ptrid_block4_lus_dgemm


! note: a(:,3:4,3,:) = 0, a(:,3:4,5,:) = 0, and a(:,1:2,4,:) = 0
! a(:,1:2,5,:) packed in a(:,1:2,4,:)
  subroutine ptrid_block4_lus_mm( a, r, n, n1, n2 )
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(4,4,4,n), intent(in) :: a
    real(kind=c_double), dimension(4,n1,n2,n), intent(inout) :: r
    integer(c_int) :: i,j,k
    do k=2,n
      do j=1,n2
        r(:,:,j,k) = r(:,:,j,k) - matmul(a(:,:,1,k),r(:,:,j,k-1))
      end do
    end do
    do k=1,n-2
      do j=1,n2
        r(:,:,j,n) = r(:,:,j,n) - matmul(a(:,1:2,4,k),r(1:2,:,j,k))
      end do
    end do
    do j=1,n2
      r(:,:,j,n) = matmul(a(:,:,2,n),r(:,:,j,n))
    end do
    do k=n-1,1,-1
      do j=1,n2
        r(:,:,j,k) = matmul( a(:,:,2,k), r(:,:,j,k)-matmul(a(:,1:2,3,k),r(1:2,:,j,k+1))-matmul(a(:,3:4,4,k),r(3:4,:,j,n)) )
      end do
    end do
  end subroutine ptrid_block4_lus_mm


  subroutine bpentLUS3x(c, r, n, n1, n2, n3)
    implicit none
    integer(c_int), intent(in) :: n,n1,n2,n3
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2,n3), intent(inout) :: r
    integer(c_int) :: i,j,k
    if( n > n1 ) return
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do k=1,n3
     do j=1,n2
      do i=1,n-2
        r(i+1,j,k) = r(i+1,j,k) - r(i,j,k)*c(i+1,2)
        r(i+2,j,k) = r(i+2,j,k) - r(i,j,k)*c(i+2,1)
      end do
      r(n,j,k) = (r(n,j,k) - r(n-1,j,k)*c(n,2))*c(n,3)
      r(n-1,j,k) = (r(n-1,j,k) - c(n-1,4)*r(n,j,k))*c(n-1,3)
      do i=n-2,1,-1
        r(i,j,k) = (r(i,j,k) - c(i,4)*r(i+1,j,k) - c(i,5)*r(i+2,j,k))*c(i,3)
      end do
     end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine bpentLUS3x


  SUBROUTINE ppentLUS3x(c,r,n,n1,n2,n3)
    IMPLICIT NONE
    integer(c_int), intent(in) :: n,n1,n2,n3
    real(kind=c_double), DIMENSION(n,9), INTENT(IN) :: c
    real(kind=c_double), DIMENSION(n1,n2,n3), INTENT(INOUT) :: r
    real(kind=c_double) :: tmp1,tmp2
    integer(c_int) :: i,j,k
    if( n > n1 ) return
    !$omp target teams distribute parallel do collapse(2) private(tmp1,tmp2) if(gpu_kernel==1) !$fexl-async.sync_var. 
    do k=1,n3
     do j=1,n2
      r(2,j,k)=r(2,j,k)-c(2,2)*r(1,j,k)
      tmp1 = c(1,6)*r(1,j,k)+c(2,6)*r(2,j,k)
      tmp2 = c(1,7)*r(1,j,k)+c(2,7)*r(2,j,k)
      DO i=3,N-2
        r(i,j,k)=r(i,j,k)-(c(i,2)*r(i-1,j,k)+c(i,1)*r(i-2,j,k))
        tmp1 = tmp1+c(i,6)*r(i,j,k)
        tmp2 = tmp2+c(i,7)*r(i,j,k)
      END DO
      r(N-1,j,k)=r(N-1,j,k)-tmp1
      r(N,j,k)=(r(N,j,k)-tmp2-c(N-1,7)*r(N-1,j,k))*c(N,3)
      r(N-1,j,k)=(r(N-1,j,k)-c(N-1,9)*r(N,j,k))*c(N-1,3)
      r(N-2,j,k)=(r(N-2,j,k)-c(N-2,8)*r(N-1,j,k)-c(N-2,9)*r(N,j,k))*c(N-2,3)
      r(N-3,j,k)=(r(N-3,j,k)-(c(N-3,4)*r(N-2,j,k)+c(N-3,8)*r(N-1,j,k)+c(N-3,9)*r(N,j,k)))*c(N-3,3)
      DO i=N-4,1,-1
        r(i,j,k)=(r(i,j,k)-(c(i,4)*r(i+1,j,k)+c(i,5)*r(i+2,j,k)+c(i,8)*r(n-1,j,k)+c(i,9)*r(n,j,k)))*c(i,3)
      END DO
     end do
    end do
    !$omp end target teams distribute parallel do
  END SUBROUTINE ppentLUS3x


  subroutine bpentLUS3y(c, r, n, n1, n2, n3)
    implicit none
    integer(c_int), intent(in) :: n, n1, n2, n3
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2,n3), intent(inout) :: r
    integer(c_int) :: i,j,k
    if( n > n2 ) return
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do k=1,n3
    	do i=1,n1
       	do j=1,n-2
          r(i,j+1,k) = r(i,j+1,k) - r(i,j,k)*c(j+1,2)
          r(i,j+2,k) = r(i,j+2,k) - r(i,j,k)*c(j+2,1)
        end do ! end j
        r(i,n,k) = (r(i,n,k) - r(i,n-1,k)*c(n,2))*c(n,3)
        r(i,n-1,k) = (r(i,n-1,k) - c(n-1,4)*r(i,n,k))*c(n-1,3)
      	do j=n-2,1,-1 
        	r(i,j,k) = (r(i,j,k) - c(j,4)*r(i,j+1,k) - c(j,5)*r(i,j+2,k))*c(j,3)
        end do ! end j
      end do ! end i
    end do ! end k
    !$omp end target teams distribute parallel do
  end subroutine bpentLUS3y


  subroutine bpentLUS2y(c, r, n, n1, n2)
    implicit none
    integer(c_int), intent(in) :: n, n1, n2
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2), intent(inout) :: r
    integer(c_int) :: i,j,k
    if( n > n2 ) return
      do j=1,n-2
        do i=1,n1
          r(i,j+1) = r(i,j+1) - r(i,j)*c(j+1,2)
          r(i,j+2) = r(i,j+2) - r(i,j)*c(j+2,1)
        end do
      end do
      do i=1,n1
        r(i,n) = (r(i,n) - r(i,n-1)*c(n,2))*c(n,3)
        r(i,n-1) = (r(i,n-1) - c(n-1,4)*r(i,n))*c(n-1,3)
      end do
      do j=n-2,1,-1 
       do i=1,n1
          r(i,j) = (r(i,j) - c(j,4)*r(i,j+1) - c(j,5)*r(i,j+2))*c(j,3)
        end do
      end do
  end subroutine bpentLUS2y


  SUBROUTINE ppentLUS3y(c, r, n, n1, n2, n3)
    IMPLICIT NONE
    integer(c_int), intent(in) :: n, n1, n2, n3
    real(kind=c_double), DIMENSION(n,9), INTENT(IN) :: c
    real(kind=c_double), DIMENSION(n1,n2,n3), INTENT(INOUT) :: r
    real(kind=c_double) :: tmp1,tmp2
    integer(c_int) :: i,j,k
    if( n > n2 ) return
    !$omp target teams distribute parallel do collapse(2) private(tmp1,tmp2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do k=1,n3
      do i=1,n1
        r(i,2,k) = r(i,2,k)-c(2,2)*r(i,1,k)
        tmp1 = c(1,6)*r(i,1,k)+c(2,6)*r(i,2,k)
        tmp2 = c(1,7)*r(i,1,k)+c(2,7)*r(i,2,k)
      DO j=3,N-2
          r(i,j,k) = r(i,j,k)-(c(j,2)*r(i,j-1,k)+c(j,1)*r(i,j-2,k))
          tmp1 = tmp1+c(j,6)*r(i,j,k)
          tmp2 = tmp2+c(j,7)*r(i,j,k)
      END DO
        r(i,N-1,k)=r(i,N-1,k)-tmp1
        r(i,N,k)  =(r(i,N,k)-tmp2-c(N-1,7)*r(i,N-1,k))*c(N,3)
        r(i,N-1,k)=(r(i,N-1,k)-c(N-1,9)*r(i,N,k))*c(N-1,3)
        r(i,N-2,k)=(r(i,N-2,k)-c(N-2,8)*r(i,N-1,k)-c(N-2,9)*r(i,N,k))*c(N-2,3)
        r(i,N-3,k)=(r(i,N-3,k)-(c(N-3,4)*r(i,N-2,k)+c(N-3,8)*r(i,N-1,k)+c(N-3,9)*r(i,N,k)))*c(N-3,3)
      DO j=N-4,1,-1
          r(i,j,k)=(r(i,j,k)-(c(j,4)*r(i,j+1,k)+c(j,5)*r(i,j+2,k)+c(j,8)*r(i,n-1,k)+c(j,9)*r(i,n,k)))*c(j,3)
        end do
      END DO
    end do
    !$omp end target teams distribute parallel do
  END SUBROUTINE ppentLUS3y


  subroutine bpentLUS3z(c, r, n, n1, n2, n3)
    implicit none
    integer(c_int), intent(in) :: n, n1, n2, n3
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2,n3), intent(inout) :: r
    integer(c_int) :: i,j,k
    if( n > n3 ) return
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1) !$fexl-async.sync_var.
    do j=1,n2
     do i=1,n1
    	do k=1,n-2
       r(i,j,k+1) = r(i,j,k+1) - r(i,j,k)*c(k+1,2)
       r(i,j,k+2) = r(i,j,k+2) - r(i,j,k)*c(k+2,1)
    	end do ! k
		  r(i,j,n) = (r(i,j,n) - r(i,j,n-1)*c(n,2))*c(n,3)
		  r(i,j,n-1) = (r(i,j,n-1) - c(n-1,4)*r(i,j,n))*c(n-1,3)
		  do k=n-2,1,-1
		    r(i,j,k) = (r(i,j,k) - c(k,4)*r(i,j,k+1) - c(k,5)*r(i,j,k+2))*c(k,3)
		  end do
     end do ! i
    end do ! j
    !$omp end target teams distribute parallel do
  end subroutine bpentLUS3z


  SUBROUTINE ppentLUS3z(c, r, n, n1, n2, n3)
      IMPLICIT NONE
      integer(c_int), intent(in) :: n, n1, n2, n3
      real(kind=c_double), DIMENSION(n,9), INTENT(IN) :: c
      real(kind=c_double), DIMENSION(n1,n2,n3), INTENT(INOUT) :: r
      real(kind=c_double) :: tmp1,tmp2
      integer(c_int) :: i,j,k
      if( n > n3 ) return
      !$omp target teams distribute parallel do collapse(2) private(tmp1,tmp2) if(gpu_kernel==1) !$fexl-async.sync_var. 
      do j=1,n2
       do i=1,n1
        r(i,j,2)=r(i,j,2)-c(2,2)*r(i,j,1)
        tmp1 = c(1,6)*r(i,j,1)+c(2,6)*r(i,j,2)
        tmp2 = c(1,7)*r(i,j,1)+c(2,7)*r(i,j,2)
        DO k=3,N-2
         r(i,j,k)=r(i,j,k)-(c(k,2)*r(i,j,k-1)+c(k,1)*r(i,j,k-2))
         tmp1 = tmp1+c(k,6)*r(i,j,k) 
         tmp2 = tmp2+c(k,7)*r(i,j,k) 
        END DO
        r(i,j,N-1)=r(i,j,N-1)-tmp1
        r(i,j,N)=(r(i,j,N)-tmp2-c(N-1,7)*r(i,j,N-1))*c(N,3)
        r(i,j,N-1)=(r(i,j,N-1)-c(N-1,9)*r(i,j,N))*c(N-1,3)
        r(i,j,N-2)=(r(i,j,N-2)-c(N-2,8)*r(i,j,N-1)-c(N-2,9)*r(i,j,N))*c(N-2,3)
        r(i,j,N-3)=(r(i,j,N-3)-(c(N-3,4)*r(i,j,N-2)+c(N-3,8)*r(i,j,N-1)+c(N-3,9)*r(i,j,N)))*c(N-3,3)
        DO k=N-4,1,-1
         r(i,j,k)=(r(i,j,k)-(c(k,4)*r(i,j,k+1)+c(k,5)*r(i,j,k+2)+c(k,8)*r(i,j,n-1)+c(k,9)*r(i,j,n)))*c(k,3)
        END DO ! k
       end do ! i
      end do ! j
      !$omp end target teams distribute parallel do
  END SUBROUTINE ppentLUS3z

  SUBROUTINE ppentLUS(opindex,use_ppent_opt,clu,FRHS)
! clu (input) contains the LU decomposition coefficients from ppentLUD
! FRHS contains in RHS lots on input and the solution lots on output
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: opindex
   LOGICAL(c_bool), INTENT(IN) :: use_ppent_opt
   real(kind=c_double), DIMENSION(:,:), INTENT(IN) :: clu
   real(kind=c_double), DIMENSION(:,:,:), INTENT(INOUT) :: FRHS
   INTEGER :: N,K
   INTEGER :: n1, n2, n3
!   INTEGER :: flopCount
   N = SIZE(clu,1)
   IF( SIZE(FRHS,opindex) < N ) THEN
     WRITE(*,*) 'error in ppentLUS'
     STOP
   ENDIF
   n1 = size(FRHS,1)      
   n2 = size(FRHS,2)      
   n3 = size(FRHS,3)      
   SELECT CASE(opindex)
   CASE(1) ! operate on first index
!      flopCount = n2*n3*(17*n1 - 40)
     ! First pass.  Solve L*Y = FRHS.
     ! Y(1) = FRHS(1) initially, so is already correct.
     ! Get Y(2) = FRHS(2) - Q(2)*Y(1).
      FRHS(2,:,:)=FRHS(2,:,:)-clu(2,2)*FRHS(1,:,:)
     ! Get the Y vectors for K = 3 to N-2 from Y(K) = FRHS(K) - Q(K)*Y(K-1) - P(K)*Y(K-2).
      DO K=3,N-2
        FRHS(K,:,:)=FRHS(K,:,:)-(clu(K,2)*FRHS(K-1,:,:)+clu(K,1)*FRHS(K-2,:,:))
      END DO
     ! Get Y(N-1) = FRHS(N-1) - Sum_{K=1}^{N-2} F(K)*FRHS(K).
      DO K=1,N-2
        FRHS(N-1,:,:)=FRHS(N-1,:,:)-clu(K,6)*FRHS(K,:,:)
      END DO
     ! Get Y(N) = FRHS(N) - Sum_{K=1}^{N-1} G(K)*FRHS(K).
      DO K=1,N-1
        FRHS(N,:,:)=FRHS(N,:,:)-clu(K,7)*FRHS(K,:,:)
      END DO
     ! Second pass.  Solve U*X = Y.
     ! Solve R(N)*X(N) = Y(N) for X(N).
      FRHS(N,:,:)=FRHS(N,:,:)*clu(N,3)
     ! Solve R(N-1)*X(N-1) = Y(N-1) - U(N-1)*X(N) for X(N-1).
      FRHS(N-1,:,:)=(FRHS(N-1,:,:)-clu(N-1,9)*FRHS(N,:,:))*clu(N-1,3)
     ! Solve R(N-2)*X(N-2) = Y(N-2) - T(N-2)*X(N-1) - U(N-2)*X(N) for X(N-2).
      FRHS(N-2,:,:)=(FRHS(N-2,:,:)-clu(N-2,8)*FRHS(N-1,:,:)             &
     &             -clu(N-2,9)*FRHS(N,:,:))*clu(N-2,3)
     ! Solve R(N-3)*X(N-3) = Y(N-3) - S(N-3)*X(N-2) - T(N-3)*X(N-1) - U(N-3)*X(N) for X(N-3).
      FRHS(N-3,:,:)=(FRHS(N-3,:,:)-(clu(N-3,4)*FRHS(N-2,:,:)            &
     &  +clu(N-3,8)*FRHS(N-1,:,:)+clu(N-3,9)*FRHS(N,:,:)))*clu(N-3,3)
     ! Loop over remaining rows.
     ! Solve R(K)*X(K) = Y(K) - S(K)*X(K+1) - E(K)*X(K+2) - T(K)*X(N-1) - U(K)*X(N) for X(K).
      DO K=N-4,1,-1
        FRHS(K,:,:)=(FRHS(K,:,:)-(clu(K,4)*FRHS(K+1,:,:)                &
     &    +clu(K,5)*FRHS(K+2,:,:)+clu(K,8)*FRHS(N-1,:,:)                &
     &    +clu(K,9)*FRHS(N,:,:)))*clu(K,3)
      END DO
   CASE(2) ! operate on second index
      if (use_ppent_opt) then
         call ppentlus_f77(n1,n2,n3,clu(1,1),frhs(1,1,1))
         return
      end if
!      flopCount = n1*n3*(17*n2 - 40)
     ! First pass.  Solve L*Y = FRHS.
     ! Y(1) = FRHS(1) initially, so is already correct.
     ! Get Y(2) = FRHS(2) - Q(2)*Y(1).
      FRHS(:,2,:)=FRHS(:,2,:)-clu(2,2)*FRHS(:,1,:)
     ! Get the Y vectors for K = 3 to N-2 from Y(K) = FRHS(K) - Q(K)*Y(K-1) - P(K)*Y(K-2).
      DO K=3,N-2
        FRHS(:,K,:)=FRHS(:,K,:)-(clu(K,2)*FRHS(:,K-1,:)+clu(K,1)*FRHS(:,K-2,:))
      END DO
     ! Get Y(N-1) = FRHS(N-1) - Sum_{K=1}^{N-2} F(K)*FRHS(K).
      DO K=1,N-2
        FRHS(:,N-1,:)=FRHS(:,N-1,:)-clu(K,6)*FRHS(:,K,:)
      END DO
     ! Get Y(N) = FRHS(N) - Sum_{K=1}^{N-1} G(K)*FRHS(K).
      DO K=1,N-1
        FRHS(:,N,:)=FRHS(:,N,:)-clu(K,7)*FRHS(:,K,:)
      END DO
     ! Second pass.  Solve U*X = Y.
     ! Solve R(N)*X(N) = Y(N) for X(N).
      FRHS(:,N,:)=FRHS(:,N,:)*clu(N,3)
     ! Solve R(N-1)*X(N-1) = Y(N-1) - U(N-1)*X(N) for X(N-1).
      FRHS(:,N-1,:)=(FRHS(:,N-1,:)-clu(N-1,9)*FRHS(:,N,:))*clu(N-1,3)
     ! Solve R(N-2)*X(N-2) = Y(N-2) - T(N-2)*X(N-1) - U(N-2)*X(N) for X(N-2).
      FRHS(:,N-2,:)=(FRHS(:,N-2,:)-clu(N-2,8)*FRHS(:,N-1,:)             &
     &             -clu(N-2,9)*FRHS(:,N,:))*clu(N-2,3)
     ! Solve R(N-3)*X(N-3) = Y(N-3) - S(N-3)*X(N-2) - T(N-3)*X(N-1) - U(N-3)*X(N) for X(N-3).
      FRHS(:,N-3,:)=(FRHS(:,N-3,:)-(clu(N-3,4)*FRHS(:,N-2,:)            &
     &  +clu(N-3,8)*FRHS(:,N-1,:)+clu(N-3,9)*FRHS(:,N,:)))*clu(N-3,3)
     ! Loop over remaining rows.
     ! Solve R(K)*X(K) = Y(K) - S(K)*X(K+1) - E(K)*X(K+2) - T(K)*X(N-1) - U(K)*X(N) for X(K).
      DO K=N-4,1,-1
        FRHS(:,K,:)=(FRHS(:,K,:)-(clu(K,4)*FRHS(:,K+1,:)                &
     &    +clu(K,5)*FRHS(:,K+2,:)+clu(K,8)*FRHS(:,N-1,:)                &
     &    +clu(K,9)*FRHS(:,N,:)))*clu(K,3)
      END DO
   CASE(3) ! operate on third index
      if (use_ppent_opt) then
         call ppentlus_f77(n1*n2,n3,1,clu(1,1),frhs(1,1,1))
         return
      end if
!      flopCount = n1*n2*(17*n3 - 40)
     ! First pass.  Solve L*Y = FRHS.
     ! Y(1) = FRHS(1) initially, so is already correct.
     ! Get Y(2) = FRHS(2) - Q(2)*Y(1).
      FRHS(:,:,2)=FRHS(:,:,2)-clu(2,2)*FRHS(:,:,1)
     ! Get the Y vectors for K = 3 to N-2 from Y(K) = FRHS(K) - Q(K)*Y(K-1) - P(K)*Y(K-2).
      DO K=3,N-2
        FRHS(:,:,K)=FRHS(:,:,K)-(clu(K,2)*FRHS(:,:,K-1)+clu(K,1)*FRHS(:,:,K-2))
      END DO
     ! Get Y(N-1) = FRHS(N-1) - Sum_{K=1}^{N-2} F(K)*FRHS(K).
      DO K=1,N-2
        FRHS(:,:,N-1)=FRHS(:,:,N-1)-clu(K,6)*FRHS(:,:,K)
      END DO
     ! Get Y(N) = FRHS(N) - Sum_{K=1}^{N-1} G(K)*FRHS(K).
      DO K=1,N-1
        FRHS(:,:,N)=FRHS(:,:,N)-clu(K,7)*FRHS(:,:,K)
      END DO
     ! Second pass.  Solve U*X = Y.
     ! Solve R(N)*X(N) = Y(N) for X(N).
      FRHS(:,:,N)=FRHS(:,:,N)*clu(N,3)
     ! Solve R(N-1)*X(N-1) = Y(N-1) - U(N-1)*X(N) for X(N-1).
      FRHS(:,:,N-1)=(FRHS(:,:,N-1)-clu(N-1,9)*FRHS(:,:,N))*clu(N-1,3)
     ! Solve R(N-2)*X(N-2) = Y(N-2) - T(N-2)*X(N-1) - U(N-2)*X(N) for X(N-2).
      FRHS(:,:,N-2)=(FRHS(:,:,N-2)-clu(N-2,8)*FRHS(:,:,N-1)             &
     &             -clu(N-2,9)*FRHS(:,:,N))*clu(N-2,3)
     ! Solve R(N-3)*X(N-3) = Y(N-3) - S(N-3)*X(N-2) - T(N-3)*X(N-1) - U(N-3)*X(N) for X(N-3).
      FRHS(:,:,N-3)=(FRHS(:,:,N-3)-(clu(N-3,4)*FRHS(:,:,N-2)            &
     &  +clu(N-3,8)*FRHS(:,:,N-1)+clu(N-3,9)*FRHS(:,:,N)))*clu(N-3,3)
     ! Loop over remaining rows.
     ! Solve R(K)*X(K) = Y(K) - S(K)*X(K+1) - E(K)*X(K+2) - T(K)*X(N-1) - U(K)*X(N) for X(K).
      DO K=N-4,1,-1
        FRHS(:,:,K)=(FRHS(:,:,K)-(clu(K,4)*FRHS(:,:,K+1)                &
     &    +clu(K,5)*FRHS(:,:,K+2)+clu(K,8)*FRHS(:,:,N-1)                &
     &    +clu(K,9)*FRHS(:,:,N)))*clu(K,3)
      END DO
   END SELECT
  END SUBROUTINE ppentLUS

 END MODULE LES_pentadiagonal
 

! ===================================================================
! MLW: 9/6/2005
! This routine is an optimized version of ppentLUS for the IBM BGL.
! In reality, it probably would improve performance on other IBM
! architectures as well.  The main goals are:
! (1) Try to keep streaming through memory with stride 1 so as to
!     keep the prefetch engines feeding the L1 cache.  
! (2) Unroll the inner loops (first dimension).  The solve is acting
!     along the second dimension so the first and third dimensions are
!     independent.  Unrolling allows multiple independent floating
!     point operations to be available to the FPUs so that the pipelines
!     remain full.
! It is assumed opindex==2, but this routine can be called for
! the case where opindex is 3.  In that case, fold the first two
! dimensions of frhs into n1 and set n3=1.
!
! IMPORTANT: This routine should not have a F90 interface defined
! for it.  It should look like an external F77 routine at each point
! of call.  This allows us to play the trick of folding the first
! and second dimensions like in good old F77.  In fact, from the
! caller, clu and frhs should look like scalar arguments.  
!
! In the case of opindex=2, you would call it as follows:
!    call ppentlus_f77(n1, n2, n3, clu(1,1), frhs(1,1,1), ...)
!
! In the case of opindex=3, you would call it as follows:
!    call ppentlus_f77(n1*n2, n3, 1, clu(1,1), frhs(1,1,1), ...)
!
! It may not be the case that all compilers will allow this trick
! but the IBM compilers do.  
! Also, it is the users responsibility to KNOW that the storage of
! the calling arrays is contiguous in normal FORTRAN77 order.
! Not to worry, the code will crash and burn if these conditions
! do not hold.
! ===================================================================
 subroutine ppentlus_f77(n1,n2,n3,clu,frhs)
  USE LES_nrutil
  IMPLICIT NONE
  INTEGER   n1,n2,n3
  DOUBLE PRECISION clu(n2,9)
  DOUBLE PRECISION FRHS(n1,n2,n3)
  DOUBLE PRECISION, dimension(n1) :: sum1, sum2
!  INTEGER  flopCount
  INTEGER  i, j, k, N
!  flopCount = n1*n3*(17*n2 - 40)
  N = n2
  do k = 1, n3
     do i = 1, n1
        FRHS(i,2,k)=FRHS(i,2,k)-clu(2,2)*FRHS(i,1,k)
        sum1(i) = clu(1,6)*FRHS(i,1,k) + clu(2,6)*FRHS(i,2,k)
        sum2(i) = clu(1,7)*FRHS(i,1,k) + clu(2,7)*FRHS(i,2,k)
     end do
     ! Get the Y vectors for K = 3 to N-2 from Y(K) = FRHS(K) - Q(j)*Y(j-1) - P(j)*Y(j-2).
     DO j=3,N-2
        do i = 1, n1
           FRHS(i,j,k)=FRHS(i,j,k)-(clu(j,2)*FRHS(i,j-1,k)       &
                &      +clu(j,1)*FRHS(i,j-2,k))
           sum1(i) = sum1(i) + clu(j,6)*FRHS(i,j,k)
           sum2(i) = sum2(i) + clu(j,7)*FRHS(i,j,k)
        end do
     END DO
     do i = 1, n1
        FRHS(i,N-1,k) = FRHS(i,N-1,k) - sum1(i)
        FRHS(i,N,k) = FRHS(i,N,k) - sum2(i) - clu(N-1,7)*FRHS(i,N-1,k)
     end do
  END DO
  ! Second pass.  Solve U*X = Y.
  do k = 1, n3
     do i = 1, n1
        FRHS(i,N,k)=FRHS(i,N,k)*clu(N,3)
        ! Solve R(N-1)*X(N-1) = Y(N-1) - U(N-1)*X(N) for X(N-1).
        FRHS(i,N-1,k)=(FRHS(i,N-1,k)-clu(N-1,9)*FRHS(i,N,k))*clu(N-1,3)
        ! Solve R(N-2)*X(N-2) = Y(N-2) - T(N-2)*X(N-1) - U(N-2)*X(N) for X(N-2).
        FRHS(i,N-2,k)=(FRHS(i,N-2,k)-clu(N-2,8)*FRHS(i,N-1,k)             &
             &             -clu(N-2,9)*FRHS(i,N,k))*clu(N-2,3)
        ! Solve R(N-3)*X(N-3) = Y(N-3) - S(N-3)*X(N-2) - T(N-3)*X(N-1) - U(N-3)*X(N) for X(N-3).
        FRHS(i,N-3,k)=(FRHS(i,N-3,k)-(clu(N-3,4)*FRHS(i,N-2,k)            &
             &  +clu(N-3,8)*FRHS(i,N-1,k)+clu(N-3,9)*FRHS(i,N,k)))*clu(N-3,3)
     end do
     ! Loop over remaining rows.
     ! Solve R(j)*X(j) = Y(j) - S(j)*X(j+1) - E(j)*X(j+2) - T(j)*X(N-1) - U(j)*X(N) for X(j).
     DO j=N-4,1,-1
        do i = 1, n1
           FRHS(i,j,k)=(FRHS(i,j,k)-(clu(j,4)*FRHS(i,j+1,k)             &
                &    +clu(j,5)*FRHS(i,j+2,k)+clu(j,8)*FRHS(i,N-1,k)     &
                &    +clu(j,9)*FRHS(i,N,k)))*clu(j,3)
        end do
     END DO
  end DO
 end subroutine ppentlus_f77

