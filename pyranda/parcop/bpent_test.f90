#define SIZE 64
module pent
  use iso_c_binding
  contains
  subroutine bpentLUS3x(c, r, n, n1, n2, n3)
    implicit none
    integer(c_int), intent(in) :: n,n1,n2,n3
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2,n3), intent(inout) :: r
    integer(c_int) :: i,j,k
    if( n > n1 ) return
    do k=1,n3
    do j=1,n2
      do i=1,n-2
        r(i+1,j,k) = r(i+1,j,k) - r(i,j,k)*c(i+1,2)
        r(i+2,j,k) = r(i+2,j,k) - r(i,j,k)*c(i+2,1)
      end do
      r(n,j,k) = r(n,j,k) - r(n-1,j,k)*c(n,2)
      r(n,j,k) = r(n,j,k)*c(n,3)
      r(n-1,j,k) = (r(n-1,j,k) - c(n-1,4)*r(n,j,k))*c(n-1,3)
      do i=n-2,1,-1
        r(i,j,k) = (r(i,j,k) - c(i,4)*r(i+1,j,k) - c(i,5)*r(i+2,j,k))*c(i,3)
      end do
    end do
    end do
  end subroutine bpentLUS3x
  subroutine bpentLUS3x_alt(c, r, n, n1, n2, n3)
    implicit none
    integer(c_int), intent(in) :: n,n1,n2,n3
    real(kind=c_double), dimension(n,5), intent(in) :: c
    real(kind=c_double), dimension(n1,n2,n3), intent(inout) :: r
    integer(c_int) :: i,j,k,z
    real(kind=c_double) :: qm1, qm2, thisq !q_(i-1) and q_(i-2)

    real(kind=c_double), dimension(n,n) :: q,p,x
    real(kind=c_double), dimension(n1,n2,n3) :: rout
    if( n > n1 ) return

200 FORMAT(' ', 4F10.4)
    PRINT *, "c(1:4,1:2): "
    WRITE(*,200) c(1:4,1)
    WRITE(*,200) c(1:4,2)
    !**** Build the Q matrix ****
    !**   q(i,j) is the contribution of r(j) to the forward propagated r(i)
    !**   It can be conceptualized by drawing the contributions thusly
    !**   
    !**          c8      c7      c6      c5      c4      c3      c2
    !**          <--     <--     <--     <--     <--     <--     <--
    !**      ___/   \___/   \___/   \___/   \___/   \___/   \___/   \___  
    !**     |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   | 
    !**     | 8 |   | 7 |   | 6 |   | 5 |   | 4 |   | 3 |   | 2 |   | 1 |
    !**     |___|   |___|   |___|   |___|   |___|   |___|   |___|   |___|
    !**         \       \   /   \   /   \   /   \   /   \   /       /    
    !**          \       <-/-----\--     <-/-----\--     <-/--------
    !**           <--------  d7   <--------   d5  <--------    d3 
    !**               d8              d6              d4
    !**
    !**     So, the coefficient of the contribution of the original r2 to the new r5
    !**     is c3*c4*c5 + d4*c5 + c3*d5
    !**     If q(i,j) is the contribution of the original r(j) to the new r(i), then
    !**     there is a recurrence relation. q(i,j) = c(i)*q(i-1,j) + d(i)*q(i-2,j)
    do i=1,n
      qm1 = 1.0
      qm2 = 0.0
      do j=1,n
        if (j.lt.i) THEN
           q(i,j) = 0.0d0
        elseif (j.eq.i) THEN
           q(i,j)= 1.0d0
        else
           q(i,j) = -c(j,2)*qm1 - c(j,1)*qm2
          ! PRINT *, "q(", i, ", ",j, ") = -",c(j,2), " * ", qm1, " + ", c(j,1), " * ", qm2, " = ", q(i,j)
           qm2 = qm1
           qm1 = q(i,j)
        endif 
      enddo
    enddo

    PRINT *, "q(", SIZE-3, ":", SIZE,",1:4,1)"
    do i=1,4
      WRITE(*,200) q(SIZE-3:SIZE,i)
    enddo

    do i=1,n
      qm1 = c(i,3)
      qm2 = 0.0
      do j=1,n
        if (j.lt.i) THEN
           p(j,i) = 0.0d0
        elseif (j.eq.i) THEN
           p(j,i)= c(i,3)
        else
           p(j,i) = -c(j,3)*c(j-1,4)*qm1 - c(j,3)*c(j-2,5)*qm2
          ! PRINT *, "q(", i, ", ",j, ") = -",c(j,2), " * ", qm1, " + ", c(j,1), " * ", qm2, " = ", q(i,j)
           qm2 = qm1
           qm1 = p(j,i)
        endif 
      enddo
    enddo
    
    do i=1,n
    do j=1,n
       x(i,j) = 0
       do k=1,n
          x(i,j) = x(i,j) + q(i,k)*p(k,j)
       enddo
    enddo
    enddo
    
    !x = p

    do k=1,n3
    do j=1,n2
    do i=1,n
      thisq=0
      do z=1,n
         thisq = thisq + x(z,i)*r(z,j,k)
      enddo
      rout(i,j,k) = thisq
    enddo
    enddo
    enddo
    r = rout
  end subroutine bpentLUS3x_alt
end module pent

program main

  use iso_c_binding
  use pent
  implicit none
  
  integer(c_int) :: n1, n2, n3, i, j, k
  real(kind=c_double) :: c(SIZE,5), r(SIZE,SIZE,SIZE), ralt(SIZE,SIZE,SIZE)
  n1 = SIZE
  n2 = SIZE
  n3 = SIZE
  do i=1,n1 
  do j=1,n2 
  do k=1,n3 
    !r(i,j,k) = 1.0d0 + 0.744d0*MOD(i+j+k,37) - 0.412d0*MOD(i*j+k,7) - 0.330d0*MOD(i+j*k,19) 
    r(i,j,k) = 0.0d0
    if (i.eq.j.and.k.eq.1) THEN
       r(i,j,k) = 1.0d0
    endif
  enddo
  enddo
  enddo

  do i=1,5
  do j=1,n1
    c(j,i) = 1.0 + 0.0489d0*MOD(i*j,3) + 0.0721d0*MOD(i*i-j,7) -0.12d0*MOD(j/i+j,11)
    !c(j,i) = j*1.0d0 + 10000.0d0*i
  enddo
  enddo
  ralt = r
  !PRINT *,c
  call bpentLUS3x_alt(c, ralt, n1, n1, n2, n3)
  call bpentLUS3x(c, r, n1, n1, n2, n3)
200 FORMAT(' ', 4F10.4)
  PRINT *, "r(", SIZE-3, ":", SIZE,",1:4,1)"
  do i=1,4
    WRITE(*,200) r(SIZE-3:SIZE,i,1)
  enddo
  PRINT *, "ralt(", SIZE-3, ":", SIZE,",1:4,1)"
  do i=1,4
    WRITE(*,200) ralt(SIZE-3:SIZE,i,1)
  enddo
  
  do i=1,SIZE
  do j=1,SIZE
  do k=1,SIZE
    IF (abs(r(i,j,k) - ralt(i,j,k)).gt.0.000001) THEN 
       PRINT *, "(", i,",",j,",",k,"): ", r(i,j,k), " != ", ralt(i,j,k)
    ENDIF
  enddo
  enddo
  enddo

  IF (.false.) THEN
  PRINT *, r(6,49,18)
  PRINT *, ralt(6,49,18)
  PRINT *, r(27,9,6)
  PRINT *, ralt(27,9,6)
  PRINT *, r(1,1,1)
  PRINT *, ralt(1,1,1)
  PRINT *, r(2,1,1)
  PRINT *, ralt(2,1,1)
  PRINT *, r(1,2,1)
  PRINT *, ralt(1,2,1)
  PRINT *, r(2,2,1)
  PRINT *, ralt(2,2,1)
  PRINT *, r(1,1,2)
  PRINT *, ralt(1,1,2)
  PRINT *, r(2,1,2)
  PRINT *, ralt(2,1,2)
  PRINT *, r(1,2,2)
  PRINT *, ralt(1,2,2)
  PRINT *, r(2,2,2)
  PRINT *, ralt(2,2,2)
  PRINT *, r(58,54,12)
  PRINT *, ralt(58,54,12)
  PRINT *, r(31,10,61)
  PRINT *, ralt(31,10,61)
  ENDIF

end program
