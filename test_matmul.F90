program main
  real*8, dimension(2,3) :: a
  real*8, dimension(3,4) :: b
  integer :: i,j,k
  do i=1,2
     do j=1,3
        a(i,j) = 1
     enddo
  enddo
  do i=1,3
     do j=1,4
        b(i,j) = 2
     enddo
  enddo
  print *, matmul(a,b)
end program main
