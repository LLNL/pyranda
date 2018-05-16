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
! BLOCK_MATH MODULE
! defines simple arithmetic for 2x2 blocks and 2 vectors
! defines simple arithmetic for 4x4 blocks and 4 vectors

! the order in multiplications is important, 
!   (a*b == AB) /= (b*a == BA)

! divisions are only done from one direction currently:
!   a/b = (.inv.B)A
! not A(.inv.B) as one would think!
!===================================================================================================
 MODULE LES_blockmath
!===================================================================================================  
  integer, parameter :: block2_size=2
  integer, parameter :: block4_size=4

  type block2
	real(kind=8), dimension(block2_size,block2_size) :: array
  end type block2

  type vector2
	real(kind=8), dimension(block2_size) :: array
  end type vector2

  type block4
	real(kind=8), dimension(block4_size,block4_size) :: array
  end type block4

  type vector4
	real(kind=8), dimension(block4_size) :: array
  end type vector4

  interface assignment ( = )
	module procedure block2_gets_real, block2_gets_array, &
		block2_gets_vector2, vector2_gets_real, vector2_gets_array, &
		block4_gets_real, block4_gets_vector4, block4_gets_array, &
		vector4_gets_real, vector4_gets_array
  end interface

  interface operator ( + )
	module procedure block2_plus_block2, block2_plus_real, real_plus_block2, & 
		vector2_plus_vector2, real_plus_vector2, vector2_plus_real, &
		block4_plus_block4, block4_plus_real, real_plus_block4, & 
		vector4_plus_vector4, real_plus_vector4, vector4_plus_real
  end interface

  interface operator ( - )
	module procedure block2_minus_block2, minus_block2, block2_minus_real, real_minus_block2, & 
		vector2_minus_vector2, minus_vector2, real_minus_vector2, vector2_minus_real, &
		block4_minus_block4, minus_block4, block4_minus_real, real_minus_block4, & 
		vector4_minus_vector4, minus_vector4, real_minus_vector4, vector4_minus_real
  end interface

  interface operator ( * )
	module procedure block2_times_block2, vector2_times_vector2, block2_times_real, real_times_block2, & 
		vector2_times_real, real_times_vector2, block2_times_vector2, vector2_times_block2, &
		block4_times_block4, vector4_times_vector4, block4_times_real, real_times_block4, & 
		vector4_times_real, real_times_vector4, block4_times_vector4, vector4_times_block4
  end interface

  interface operator ( .inv. )
	module procedure block2_invert, block4_invert
  end interface

  interface operator ( / )
	module procedure block2_div_block2, block2_div_real, real_div_block2, vector2_div_block2, &
		block4_div_block4, block4_div_real, real_div_block4, vector4_div_block4
  end interface

  interface transpose
	module procedure block2_transpose, block4_transpose
  end interface
!===================================================================================================
  CONTAINS
!===================================================================================================
! Assignment

  subroutine block2_gets_real(b,r)
	type(block2), intent(out) :: b
	real(kind=8), intent(in) :: r
	b%array(1,1) = r
	b%array(2,1) = 0.0_8
	b%array(1,2) = 0.0_8
	b%array(2,2) = r
	return
  end subroutine block2_gets_real

  subroutine block2_gets_array(b,a)
	type(block2), intent(out) :: b
	real(kind=8), dimension(block2_size,block2_size), intent(in) :: a
	b%array = a
	return
  end subroutine block2_gets_array

  subroutine block2_gets_vector2(b,v)
! this spreads the vector across columns
! use transpose to get spread in the other dimension
	type(block2), intent(out) :: b
	type(vector2), intent(in) :: v
	b%array = spread( v%array, 2, block2_size )
	return
  end subroutine block2_gets_vector2

  subroutine vector2_gets_real(v,r)
	type(vector2), intent(out) :: v
	real(kind=8), intent(in) :: r
	v%array = r
	return
  end subroutine vector2_gets_real

  subroutine vector2_gets_array(v,a)
	type(vector2), intent(out) :: v
	real(kind=8), dimension(block2_size), intent(in) :: a
	v%array = a
	return
  end subroutine vector2_gets_array

  subroutine block4_gets_real(b,r)
	type(block4), intent(out) :: b
	real(kind=8), intent(in) :: r
	integer :: i
	b%array = 0.0_8
	do i=1,block4_size
	  b%array(i,i) = r
    end do
	return
  end subroutine block4_gets_real

  subroutine block4_gets_array(b,a)
	type(block4), intent(out) :: b
	real(kind=8), dimension(block4_size,block4_size), intent(in) :: a
	b%array = a
	return
  end subroutine block4_gets_array

  subroutine block4_gets_vector4(b,v)
! this spreads the vector across columns
! use transpose to get spread in the other dimension
	type(block4), intent(out) :: b
	type(vector4), intent(in) :: v
	b%array = spread( v%array, 2, block4_size )
	return
  end subroutine block4_gets_vector4

  subroutine vector4_gets_real(v,r)
	type(vector4), intent(out) :: v
	real(kind=8), intent(in) :: r
	v%array = r
	return
  end subroutine vector4_gets_real

  subroutine vector4_gets_array(v,a)
	type(vector4), intent(out) :: v
	real(kind=8), dimension(block4_size), intent(in) :: a
	v%array = a
	return
  end subroutine vector4_gets_array

! Addition

  function block2_plus_block2( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	type(block2), intent(in) :: b
	c%array = a%array+b%array
	return
  end function block2_plus_block2

  function block2_plus_real( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array(1,1) = a%array(1,1)+b
	c%array(2,1) = a%array(2,1)
	c%array(1,2) = a%array(1,2)
	c%array(2,2) = a%array(2,2)+b
	return
  end function block2_plus_real

  function real_plus_block2( a,b ) result (c)
	type(block2) :: c
	real(kind=8), intent(in) :: a
	type(block2), intent(in) :: b
	c%array(1,1) = a+b%array(1,1)
	c%array(2,1) = b%array(2,1)
	c%array(1,2) = b%array(1,2)
	c%array(2,2) = a+b%array(2,2)
	return
  end function real_plus_block2

  function vector2_plus_vector2( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = a%array+b%array
	return
  end function vector2_plus_vector2

  function vector2_plus_real( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array+b
	return
  end function vector2_plus_real

  function real_plus_vector2( a,b ) result (c)
	type(vector2) :: c
	real(kind=8), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = a+b%array
	return
  end function real_plus_vector2

  function block4_plus_block4( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	type(block4), intent(in) :: b
	c%array = a%array+b%array
	return
  end function block4_plus_block4

  function block4_plus_real( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	real(kind=8), intent(in) :: b
	integer :: i
	c%array = a%array
	do i=1,block4_size
	  c%array(i,i) = c%array(i,i)+b
    end do
	return
  end function block4_plus_real

  function real_plus_block4( a,b ) result (c)
	type(block4) :: c
	real(kind=8), intent(in) :: a
	type(block4), intent(in) :: b
	integer :: i
	c%array = b%array
	do i=1,block4_size
	  c%array(i,i) = a+c%array(i,i)
    end do
	return
  end function real_plus_block4

  function vector4_plus_vector4( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = a%array+b%array
	return
  end function vector4_plus_vector4

  function vector4_plus_real( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array+b
	return
  end function vector4_plus_real

  function real_plus_vector4( a,b ) result (c)
	type(vector4) :: c
	real(kind=8), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = a+b%array
	return
  end function real_plus_vector4

! Subtraction

  function block2_minus_block2( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	type(block2), intent(in) :: b
	c%array = a%array-b%array
	return
  end function block2_minus_block2

  function block2_minus_real( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array(1,1) = a%array(1,1)-b
	c%array(2,1) = a%array(2,1)
	c%array(1,2) = a%array(1,2)
	c%array(2,2) = a%array(2,2)-b
	return
  end function block2_minus_real

  function real_minus_block2( a,b ) result (c)
	type(block2) :: c
	real(kind=8), intent(in) :: a
	type(block2), intent(in) :: b
	c%array(1,1) = a-b%array(1,1)
	c%array(2,1) = -b%array(2,1)
	c%array(1,2) = -b%array(1,2)
	c%array(2,2) = a-b%array(2,2)
	return
  end function real_minus_block2

  function vector2_minus_vector2( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = a%array-b%array
	return
  end function vector2_minus_vector2

  function vector2_minus_real( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array-b
	return
  end function vector2_minus_real

  function real_minus_vector2( a,b ) result (c)
	type(vector2) :: c
	real(kind=8), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = a-b%array
	return
  end function real_minus_vector2

  function minus_vector2( a ) result (b)
	type(vector2) :: b
	type(vector2), intent(in) :: a
	b%array = -a%array
	return
  end function minus_vector2

  function minus_block2( a ) result (b)
	type(block2) :: b
	type(block2), intent(in) :: a
	b%array = -a%array
	return
  end function minus_block2

  function block4_minus_block4( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	type(block4), intent(in) :: b
	c%array = a%array-b%array
	return
  end function block4_minus_block4

  function block4_minus_real( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	real(kind=8), intent(in) :: b
	integer :: i
	c%array = a%array
	do i=1,block4_size
	  c%array(i,i) = c%array(i,i)-b
    end do
	return
  end function block4_minus_real

  function real_minus_block4( a,b ) result (c)
	type(block4) :: c
	real(kind=8), intent(in) :: a
	type(block4), intent(in) :: b
	integer :: i
	c%array = b%array
	do i=1,block4_size
	  c%array(i,i) = a-c%array(i,i)
    end do
	return
  end function real_minus_block4

  function vector4_minus_vector4( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = a%array-b%array
	return
  end function vector4_minus_vector4

  function vector4_minus_real( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array-b
	return
  end function vector4_minus_real

  function real_minus_vector4( a,b ) result (c)
	type(vector4) :: c
	real(kind=8), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = a-b%array
	return
  end function real_minus_vector4

  function minus_vector4( a ) result (b)
	type(vector4) :: b
	type(vector4), intent(in) :: a
	b%array = -a%array
	return
  end function minus_vector4

  function minus_block4( a ) result (b)
	type(block4) :: b
	type(block4), intent(in) :: a
	b%array = -a%array
	return
  end function minus_block4

! Multiplication

  function block2_times_block2( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	type(block2), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function block2_times_block2

  function vector2_times_vector2( a,b ) result (c)
	real(kind=8) :: c
	type(vector2), intent(in) :: a
	type(vector2), intent(in) :: b
	c = dot_product(a%array,b%array)
	return
  end function vector2_times_vector2

  function block2_times_vector2( a,b ) result (c)
	type(vector2) :: c
	type(block2), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function block2_times_vector2

  function vector2_times_block2( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	type(block2), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function vector2_times_block2

  function block2_times_real( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = b*a%array
	return
  end function block2_times_real

  function real_times_block2( a,b ) result (c)
	type(block2) :: c
	real(kind=8), intent(in) :: a
	type(block2), intent(in) :: b
	c%array = a*b%array
	return
  end function real_times_block2

  function vector2_times_real( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = b*a%array
	return
  end function vector2_times_real

  function real_times_vector2( a,b ) result (c)
	type(vector2) :: c
	real(kind=8), intent(in) :: a
	type(vector2), intent(in) :: b
	c%array = a*b%array
	return
  end function real_times_vector2

  function block4_times_block4( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	type(block4), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function block4_times_block4

  function vector4_times_vector4( a,b ) result (c)
	real(kind=8) :: c
	type(vector4), intent(in) :: a
	type(vector4), intent(in) :: b
	c = dot_product(a%array,b%array)
	return
  end function vector4_times_vector4

  function block4_times_vector4( a,b ) result (c)
	type(vector4) :: c
	type(block4), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function block4_times_vector4

  function vector4_times_block4( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	type(block4), intent(in) :: b
	c%array = matmul(a%array,b%array)
	return
  end function vector4_times_block4

  function block4_times_real( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = b*a%array
	return
  end function block4_times_real

  function real_times_block4( a,b ) result (c)
	type(block4) :: c
	real(kind=8), intent(in) :: a
	type(block4), intent(in) :: b
	c%array = a*b%array
	return
  end function real_times_block4

  function vector4_times_real( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = b*a%array
	return
  end function vector4_times_real

  function real_times_vector4( a,b ) result (c)
	type(vector4) :: c
	real(kind=8), intent(in) :: a
	type(vector4), intent(in) :: b
	c%array = a*b%array
	return
  end function real_times_vector4

! Division

  function block2_invert( a ) result (b)
	type(block2) :: b
	type(block2), intent(in) :: a
	b%array(2,2) = 1./(a%array(1,1)*a%array(2,2)-a%array(2,1)*a%array(1,2))
	b%array(1,1) = a%array(2,2)*b%array(2,2)
	b%array(2,1) = -a%array(2,1)*b%array(2,2)
	b%array(1,2) = -a%array(1,2)*b%array(2,2)
	b%array(2,2) = a%array(1,1)*b%array(2,2)
	return
  end function block2_invert

  function block2_div_block2( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	type(block2), intent(in) :: b
	c = (.inv.b) * a	! note the reverse order!
	return
  end function block2_div_block2

  function vector2_div_block2( a,b ) result (c)
	type(vector2) :: c
	type(vector2), intent(in) :: a
	type(block2), intent(in) :: b
	c = (.inv.b) * a	! note the reverse order!
	return
  end function vector2_div_block2

  function block2_div_real( a,b ) result (c)
	type(block2) :: c
	type(block2), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array / b
	return
  end function block2_div_real

  function real_div_block2( a,b ) result (c)
	type(block2) :: c
	real(kind=8), intent(in) :: a
	type(block2), intent(in) :: b
	c = (.inv.b) * a
	return
  end function real_div_block2

  function block4_invert(aa) result (bb)
    type(block4) :: bb
    type(block4), intent(in) :: aa
    type(block2) :: a,b,c,d
    a%array = aa%array(1:2,1:2)
    b%array = aa%array(1:2,3:4)
    c%array = aa%array(3:4,1:2)
    d%array = aa%array(3:4,3:4)
    d = .inv.d
    c = -d*c
    a = .inv.(a+b*c)
    b = -a*b*d
    d = d+c*b
    c = c*a
    bb%array(1:2,1:2) = a%array
    bb%array(1:2,3:4) = b%array
    bb%array(3:4,1:2) = c%array
    bb%array(3:4,3:4) = d%array
	return
  end function block4_invert

  function block4_div_block4( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	type(block4), intent(in) :: b
	c = (.inv.b) * a	! note the reverse order!
	return
  end function block4_div_block4

  function vector4_div_block4( a,b ) result (c)
	type(vector4) :: c
	type(vector4), intent(in) :: a
	type(block4), intent(in) :: b
	c = (.inv.b) * a	! note the reverse order!
	return
  end function vector4_div_block4

  function block4_div_real( a,b ) result (c)
	type(block4) :: c
	type(block4), intent(in) :: a
	real(kind=8), intent(in) :: b
	c%array = a%array / b
	return
  end function block4_div_real

  function real_div_block4( a,b ) result (c)
	type(block4) :: c
	real(kind=8), intent(in) :: a
	type(block4), intent(in) :: b
	c = (.inv.b) * a
	return
  end function real_div_block4

! Transposition

  function block2_transpose(a) result (b)
	type(block2) :: b
	type(block2), intent(in) :: a
	b%array = transpose(a%array)
	return
  end function block2_transpose

  function block4_transpose(a) result (b)
	type(block4) :: b
	type(block4), intent(in) :: a
	b%array = transpose(a%array)
	return
  end function block4_transpose
!===================================================================================================	
 END MODULE LES_blockmath
!===================================================================================================
