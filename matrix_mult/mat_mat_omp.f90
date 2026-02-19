program matrix
  use omp_lib
  implicit none

  integer :: n
  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  double precision :: t0, t1, elapsed

  n = 4096

  allocate(A(n,n), B(n,n), C(n,n))

  ! Initialize input matrices
  call random_number(B)
  call random_number(C)

  t0 = omp_get_wtime()
  call matrix_mul_omp(A, B, C, n)
  t1 = omp_get_wtime()

  elapsed = t1 - t0
  print *, 'Time elapsed is ', elapsed
  !print *, 'Checksum:', sum(A)

  deallocate(A, B, C)
end program matrix


subroutine matrix_mul_omp(A, B, C, n)
  use omp_lib
  implicit none

  integer, intent(in) :: n
  double precision, intent(in)  :: B(n,n), C(n,n)
  double precision, intent(out) :: A(n,n)

  integer :: i, j, k

  A = 0.0d0

  !$omp parallel do private (i,j,k) 
  do i = 1, n ! make this loop run in parallel over N workers
     do j = 1, n
        do k = 1, n
           A(i,j) = A(i,j) + B(i,k) * C(k,j)
        end do
     end do
  end do
  !$omp end parallel do

end subroutine matrix_mul_omp
