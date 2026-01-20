program matrix
  implicit none          ! make sure that all variables are declared
  integer :: i,j,k,n
  double precision, allocatable, dimension(:,:) :: A,B,C
  double precision :: elapsed,t0,t1
  n = 3000               ! for now, fix the matrix size
  allocate(A(n,n),B(n,n),C(n,n)) ! allocate memory space for the arrays

  ! fill the input apprays with random numbers (to avoid sparse matrices)
  call random_number(B)
  call random_number(C)
  call cpu_time(t0)        ! returns the current system clock time
  call matrix_mul(A,B,C,n) ! compute the product
  call cpu_time(t1)
  elapsed = t1-t0
  print *,'Time elapsed is ',elapsed

  deallocate(A,B,C)        ! free the memory space
end program matrix

subroutine matrix_mul(A,B,C,n)
  integer :: i,j,k
  double precision :: A(n,n),B(n,n),C(n,n)

  A = 0.0                  ! set every element of A to 0
! now do the actual work: the triple loop:
  do i=1,n
     do j=1,n
        do k=1,n
           A(i,j) = A(i,j) + B(i,k)*C(k,j)
        end do
     end do
  end do

  return
end subroutine matrix_mul
