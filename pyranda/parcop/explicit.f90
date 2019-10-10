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
module LES_explicit
  use iso_c_binding

  ! Stencil data
  
  !  First derivative, 6th order, explicit
  real(c_double), parameter :: c6e =  1.0D0 / 60.0D0
  real(c_double), parameter :: b6e = -9.0D0 / 60.0D0
  real(c_double), parameter :: a6e = 45.0D0 / 60.0D0
  real(c_double), dimension(3) :: d16e   = [a6e,b6e,c6e]
  real(c_double), dimension(5) :: d16eb3 = [1.0,-8.0,0.0,8.0,-1.0]/12.0D0
  real(c_double), dimension(4) :: d16eb2 = [-2.0,-3.0,6.0,-1.0]/12.0D0
  real(c_double), dimension(3) :: d16eb1 = [-3.0,4.0,-1.0]/2.0


  !  Fourth derivative, 4th order, explicit (for AFLES terms)
  real(c_double), parameter :: ad1e4 = 28.0D0/3.0D0
  real(c_double), parameter :: bd1e4 = -6.5D0
  real(c_double), parameter :: cd1e4 = 2.0D0
  real(c_double), parameter :: dd1e4 = -1.0D0/6.0D0
  real(c_double), dimension(4) :: d4e4 = [ad1e4,bd1e4,cd1e4,dd1e4]
  real(c_double), dimension(5) :: d4e4b3 = [1.0, -4.0, 6.0, -4.0, 1.0]
  ! b2,b1 use values of b3
  

contains

!  subroutine divVe(fx,fy,fz,df)
!    IMPLICIT NONE
!    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
!    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: df
!    DOUBLE PRECISION, DIMENSION(SIZE(df,1),SIZE(df,2),SIZE(df,3)) :: fA,fB,fC,tmp

    ! Get boundary data      
!    CALL ddx(fx,fA,patch_ptr%isymX)
!    CALL ddy(fy,fB,patch_ptr%isymY)
!    CALL ddz(fz,fC,patch_ptr%isymZ)
    ! Add up directions
!    df = fA + fB + fC       
!  end subroutine divVe

  
  subroutine der1e6(f,df,dx,direction,bc1,bc2)
    ! Assumes boundary data has been copied for periodic flow
    real(c_double), dimension(:,:,:),intent(in)  :: f
    real(c_double), dimension(:,:,:),intent(out) :: df
    real(c_double), intent(in) :: dx
    integer(c_int), intent(in) :: direction,bc1,bc2
    
    integer :: i,j,k
    double precision :: invdx
    
    invdx = 1.0 / dx
    
    
    select case(direction)
       
    case(1)

       ! interior points
       do i=4, size(f,1)-3
          do j=1,size(f,2)
             do k=1,size(f,3)
                
                df(i,j,k) = ( & 
                     d16e(3)*( f(i+3,j,k) - f(i-3,j,k) ) + &
                     d16e(2)*( f(i+2,j,k) - f(i-2,j,k) ) + &
                     d16e(1)*( f(i+1,j,k) - f(i-1,j,k) ) ) * invdx
             end do
          end do
       end do
       
       ! bcs, 0-skip
       !      1-sym
       !      2-telescoped
       !     -1-antisym
       
       if ( bc1 == 2 ) then ! standard telescoped
          do j=1,size(f,2)
             do k=1,size(f,3)
                df(1,j,k) = ( &
                     f(1,j,k)*d16eb1(1) + &
                     f(2,j,k)*d16eb1(2) + &
                     f(3,j,k)*d16eb1(3) ) * invdx
                
                df(2,j,k) = ( &
                     f(1,j,k)*d16eb2(1) + &
                     f(2,j,k)*d16eb2(2) + &
                     f(3,j,k)*d16eb2(3) + &
                     f(4,j,k)*d16eb2(4) ) * invdx
                
                df(3,j,k) = ( &
                     f(1,j,k)*d16eb3(1) + &
                     f(2,j,k)*d16eb3(2) + &
                     f(3,j,k)*d16eb3(3) + &
                     f(4,j,k)*d16eb3(4) + &
                     f(5,j,k)*d16eb3(5) ) * invdx
             end do
          end do
          
       elseif (bc1 == 1) then ! symmetry
          
          do j=1,size(f,2)
             do k=1,size(f,3)
                
                df(1,j,k) = 0.0D0
                
                df(2,j,k) = ( & 
                     d16e(3)*( f(5,j,k) - f(3,j,k) ) + &
                     d16e(2)*( f(4,j,k) - f(2,j,k) ) + &
                     d16e(1)*( f(3,j,k) - f(1,j,k) ) ) * invdx
                
                df(3,j,k) = ( & 
                     d16e(3)*( f(6,j,k) - f(2,j,k) ) + &
                     d16e(2)*( f(5,j,k) - f(1,j,k) ) + &
                     d16e(1)*( f(4,j,k) - f(2,j,k) ) ) * invdx
                
             end do
          end do
          
       elseif (bc1 == -1) then ! anti-symmtry
          
          
          do j=1,size(f,2)
             do k=1,size(f,3)
                
                df(1,j,k) = ( & 
                     d16e(3)*( f(4,j,k) + f(4,j,k) ) + &
                     d16e(2)*( f(3,j,k) + f(3,j,k) ) + &
                     d16e(1)*( f(2,j,k) + f(2,j,k) ) ) * invdx
                
                df(2,j,k) = ( & 
                     d16e(3)*( f(5,j,k) + f(3,j,k) ) + &
                     d16e(2)*( f(4,j,k) + f(2,j,k) ) + &
                     d16e(1)*( f(3,j,k) - f(1,j,k) ) ) * invdx
                
                df(3,j,k) = ( & 
                     d16e(3)*( f(6,j,k) + f(2,j,k) ) + &
                     d16e(2)*( f(5,j,k) - f(1,j,k) ) + &
                     d16e(1)*( f(4,j,k) - f(2,j,k) ) ) * invdx
                
             end do
          end do
          
       end if
                   
    case(2)

       ! interior points
       do j=4, size(f,2)-3
          do i=1,size(f,1)
             do k=1,size(f,3)
                
                df(i,j,k) = ( & 
                     d16e(3)*( f(i,j+3,k) - f(i,j-3,k) ) + &
                     d16e(2)*( f(i,j+2,k) - f(i,j-2,k) ) + &
                     d16e(1)*( f(i,j+1,k) - f(i,j-1,k) ) ) * invdx
             end do
          end do
       end do
       
       ! bcs, 0-skip
       !      1-sym
       !      2-telescoped
       !     -1-antisym
       
       if ( bc1 == 2 ) then ! standard telescoped
          do i=1,size(f,1)
             do k=1,size(f,3)
                df(i,1,k) = ( &
                     f(i,1,k)*d16eb1(1) + &
                     f(i,2,k)*d16eb1(2) + &
                     f(i,3,k)*d16eb1(3) ) * invdx
                
                df(i,2,k) = ( &
                     f(i,1,k)*d16eb2(1) + &
                     f(i,2,k)*d16eb2(2) + &
                     f(i,3,k)*d16eb2(3) + &
                     f(i,4,k)*d16eb2(4) ) * invdx
                
                df(i,3,k) = ( &
                     f(i,1,k)*d16eb3(1) + &
                     f(i,2,k)*d16eb3(2) + &
                     f(i,3,k)*d16eb3(3) + &
                     f(i,4,k)*d16eb3(4) + &
                     f(i,5,k)*d16eb3(5) ) * invdx
             end do
          end do
          
       elseif (bc1 == 1) then ! symmetry
          
          do i=1,size(f,1)
             do k=1,size(f,3)
                
                df(i,1,k) = 0.0D0
                
                df(i,2,k) = ( & 
                     d16e(3)*( f(i,5,k) - f(i,3,k) ) + &
                     d16e(2)*( f(i,4,k) - f(i,2,k) ) + &
                     d16e(1)*( f(i,3,k) - f(i,1,k) ) ) * invdx
                
                df(i,2,k) = ( & 
                     d16e(3)*( f(i,6,k) - f(i,2,k) ) + &
                     d16e(2)*( f(i,5,k) - f(i,1,k) ) + &
                     d16e(1)*( f(i,4,k) - f(i,2,k) ) ) * invdx
                
             end do
          end do
          
       elseif (bc1 == -1) then ! anti-symmtry
          
          
          do i=1,size(f,1)
             do k=1,size(f,3)
                
                df(i,1,k) = ( & 
                     d16e(3)*( f(i,4,k) + f(i,4,k) ) + &
                     d16e(2)*( f(i,3,k) + f(i,3,k) ) + &
                     d16e(1)*( f(i,2,k) + f(i,2,k) ) ) * invdx
                
                df(i,2,k) = ( & 
                     d16e(3)*( f(i,5,k) + f(i,3,k) ) + &
                     d16e(2)*( f(i,4,k) + f(i,2,k) ) + &
                     d16e(1)*( f(i,3,k) - f(i,1,k) ) ) * invdx
                
                df(i,3,k) = ( & 
                     d16e(3)*( f(i,6,k) + f(i,2,k) ) + &
                     d16e(2)*( f(i,5,k) - f(i,1,k) ) + &
                     d16e(1)*( f(i,4,k) - f(i,2,k) ) ) * invdx
                
             end do
          end do
          
       end if            
       
    case(3)

       
       ! interior points
       do k=4, size(f,3)-3
          do j=1,size(f,2)
             do i=1,size(f,1)
                
                df(i,j,k) = ( & 
                     d16e(3)*( f(i,j,k+3) - f(i,j,k-3) ) + &
                     d16e(2)*( f(i,j,k+2) - f(i,j,k-2) ) + &
                     d16e(1)*( f(i,j,k+1) - f(i,j,k-1) ) ) * invdx
             end do
          end do
       end do
       
       ! bcs, 0-skip
       !      1-sym
       !      2-telescoped
       !     -1-antisym
       
       if ( bc1 == 2 ) then ! standard telescoped
          do j=1,size(f,2)
             do i=1,size(f,1)
                df(i,j,1) = ( &
                     f(i,j,1)*d16eb1(1) + &
                     f(i,j,2)*d16eb1(2) + &
                     f(i,j,3)*d16eb1(3) ) * invdx
                
                df(i,j,2) = ( &
                     f(i,j,1)*d16eb2(1) + &
                     f(i,j,2)*d16eb2(2) + &
                     f(i,j,3)*d16eb2(3) + &
                     f(i,j,4)*d16eb2(4) ) * invdx
                
                df(i,j,3) = ( &
                     f(i,j,1)*d16eb3(1) + &
                     f(i,j,2)*d16eb3(2) + &
                     f(i,j,3)*d16eb3(3) + &
                     f(i,j,4)*d16eb3(4) + &
                     f(i,j,5)*d16eb3(5) ) * invdx
             end do
          end do
          
       elseif (bc1 == 1) then ! symmetry
          
          do j=1,size(f,2)
             do i=1,size(f,1)
                
                df(i,j,1) = 0.0D0
                
                df(i,j,2) = ( & 
                     d16e(3)*( f(i,j,5) - f(i,j,3) ) + &
                     d16e(2)*( f(i,j,4) - f(i,j,2) ) + &
                     d16e(1)*( f(i,j,3) - f(i,j,1) ) ) * invdx
                
                df(i,j,3) = ( & 
                     d16e(3)*( f(i,j,6) - f(i,j,2) ) + &
                     d16e(2)*( f(i,j,5) - f(i,j,1) ) + &
                     d16e(1)*( f(i,j,4) - f(i,j,2) ) ) * invdx
                
             end do
          end do
          
       elseif (bc1 == -1) then ! anti-symmtry
          
          
          do j=1,size(f,2)
             do i=1,size(f,1)
                
                df(i,j,1) = ( & 
                     d16e(3)*( f(i,j,4) + f(i,j,4) ) + &
                     d16e(2)*( f(i,j,3) + f(i,j,3) ) + &
                     d16e(1)*( f(i,j,2) + f(i,j,2) ) ) * invdx
                
                df(i,j,2) = ( & 
                     d16e(3)*( f(i,j,5) + f(i,j,3) ) + &
                     d16e(2)*( f(i,j,4) + f(i,j,2) ) + &
                     d16e(1)*( f(i,j,3) - f(i,j,1) ) ) * invdx
                
                df(i,j,3) = ( & 
                     d16e(3)*( f(i,j,6) + f(i,j,2) ) + &
                     d16e(2)*( f(i,j,5) - f(i,j,1) ) + &
                     d16e(1)*( f(i,j,4) - f(i,j,2) ) ) * invdx
                
             end do
          end do
          
       end if
       
    end select
    
  
  end subroutine der1e6


  subroutine der4e4(f,df,dx,direction,bc1,bc2)
    ! Assumes boundary data has been copied for periodic flow
    real(c_double), dimension(:,:,:),intent(in)  :: f
    real(c_double), dimension(:,:,:),intent(out) :: df
    real(c_double), intent(in) :: dx
    integer(c_int), intent(in) :: direction,bc1,bc2
    
    integer :: i,j,k,n
    double precision :: invdx
    
    !invdx = 1.0 / dx ! NA
    
    
    select case(direction)
       
    case(1)

       ! interior points
       do i=4, size(f,1)-3
          do j=1,size(f,2)
             do k=1,size(f,3)
                df(i,j,k) = ( d4e4(1)*f(i,j,k) + & 
                     d4e4(4)*( f(i+3,j,k) + f(i-3,j,k) ) + &
                     d4e4(3)*( f(i+2,j,k) + f(i-2,j,k) ) + &
                     d4e4(2)*( f(i+1,j,k) + f(i-1,j,k) ) )
             end do
          end do
       end do


       if ( bc1 == 2 ) then
          do j=1,size(f,2)
             do k=1,size(f,3)
                df(3,j,k) = ( &
                     f(1,j,k)*d4e4b3(1) + &
                     f(2,j,k)*d4e4b3(2) + &
                     f(3,j,k)*d4e4b3(3) + &
                     f(4,j,k)*d4e4b3(4) + &
                     f(5,j,k)*d4e4b3(5) )
                df(2,j,k) = df(3,j,k)
                df(1,j,k) = df(3,j,k)          
             end do
          end do
       endif

       if ( bc2 == 2 ) then
          n = size(f,1)
          do j=1,size(f,2)
             do k=1,size(f,3)
                df(n-2,j,k) = ( &
                     f(n,j,k)*d4e4b3(1) + &
                     f(n-1,j,k)*d4e4b3(2) + &
                     f(n-2,j,k)*d4e4b3(3) + &
                     f(n-3,j,k)*d4e4b3(4) + &
                     f(n-4,j,k)*d4e4b3(5) )
                df(n-1,j,k) = df(n-2,j,k)
                df(n,j,k)   = df(n-2,j,k)          
             end do
          end do
       endif
                           

       
    case(2)

       ! interior points
       do j=4, size(f,2)-3
          df(:,j,:) = ( d4e4(1)*f(:,j,:) + & 
               d4e4(4)*( f(:,j+3,:) + f(:,j-3,:) ) + &
               d4e4(3)*( f(:,j+2,:) + f(:,j-2,:) ) + &
               d4e4(2)*( f(:,j+1,:) + f(:,j-1,:) ) ) * invdx
       end do

    case(3)
       
       ! interior points
       do k=4, size(f,3)-3
          df(:,:,k) = ( d4e4(1)*f(:,:,k) + & 
               d4e4(4)*( f(:,:,k+3) + f(:,:,k-3) ) + &
               d4e4(3)*( f(:,:,k+2) + f(:,:,k-2) ) + &
               d4e4(2)*( f(:,:,k+1) + f(:,:,k-1) ) ) * invdx
       end do
    end select
  end subroutine der4e4

  
end module LES_explicit  
  
