program matrix
  implicit none          ! make sure that all variables are declared
  integer :: n
  double precision, allocatable, dimension(:,:) :: A,B,C
  double precision :: elapsed,t0,t1
  n = 3000               ! for now, fix the matrix size
  allocate(A(n,n),B(n,n),C(n,n)) ! allocate memory space for the arrays

  ! fill the input apprays with random numbers (to avoid sparse matrices)
  ! there is no need to fill the output array A, as it will be the output matrix
  call random_number(B)
  call random_number(C)
  call cpu_time(t0)        ! returns the current system clock time
  call matrix_mul(A,B,C,n) ! compute the product
  call cpu_time(t1)

  elapsed = t1-t0
  print *,'Time elapsed is ',elapsed

  print *, 'A(1,1)=',A(1,1)  ! print one element to make sure the compiler
                                ! does not optimize away the whole computation
  print *, 'Checksum:', sum(A)  ! print the sum of all elements as a checksum
                                ! a checksum is a singlw value calculated from a larger set of data

  deallocate(A,B,C)        ! free the memory space
end program matrix

subroutine matrix_mul(A,B,C,n)
  implicit none
  integer, intent(in) :: n
  integer :: i, j, k
  double precision, intent(out) :: A(n,n)
  double precision, intent(in) :: B(n,n), C(n,n)

  A = 0.0                  ! set every element of A to 0

! now do the actual work: the triple loop:
  do k=1,n
     do j=1,n
        do i=1,n
           A(i,j) = A(i,j) + B(i,k)*C(k,j)
        end do
     end do
  end do

  return
end subroutine matrix_mul
