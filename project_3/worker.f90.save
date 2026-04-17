subroutine worker
use globals
use auxiliary
implicit none
logical :: conv
integer :: seed,ierr,recvd_tag,tag
integer, dimension(mpi_status_size) :: status
double precision :: eig
double precision, allocatable, dimension(:,:) :: A

! Receive matrix size
call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
allocate(A(n,n))

do while(.true.)
   call mpi_recv(seed,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,ierr)
   recvd_tag = status (MPI_TAG)
   if(recvd_tag.eq.EXIT_TAG)  then            ! The exit tag was received, wind down and exit.
      print *,'Worker ',proc_num,' exiting...'
      exit
   end if
   ! The seed was received and more data are requested, so compute an eigenvalue:
   call open_random_stream(seed)
   call make_random_matrix(A,seed)
   call compute_eig(A,eig,conv)
   ! Send result
   if(conv) then
      tag=0
   else
      tag=1
   end if
   call mpi_send(eig,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
   call close_random_stream()
end do

deallocate(A)
return
end subroutine worker

subroutine make_random_matrix(A,seed)
use globals
use auxiliary
implicit none
integer :: seed,ierr,i,j
double precision :: A(n,n)

do i=1,n
   ierr=vdRNGGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,1,A(i,i),0d0,2d0)
   do j=1,i-1
      ierr=vdRNGGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,1,A(i,j),0d0,1d0)
      A(j,i)=A(i,j)
   end do
end do

return
end subroutine make_random_matrix

subroutine compute_eig(A,eig,conv)
use globals
use auxiliary
implicit none
logical :: conv
integer :: ierr,M,lwork,info,i
double precision :: A(n,n),eig,eps,mu,wr(n),wi(n),B(n,n),vl(1,n),vr(1,n),eig2
double precision, allocatable, dimension(:) :: work

lwork=4*n
allocate(work(lwork)) ! Bit of memory needed by LAPACK

eps=1d-5              ! Relative convergence criterion for power iteration (bad practice to hard-code this!)
M=int(n/2)            ! Maximal number of iterations
mu=0d0
call power(A,mu,eps,M,eig,conv)
! Check if the eigenvalue is positive (they are symmetrically distributed about 0). If not, re-compute with shift:
if(conv) then
   if(eig.lt.0) then
      mu=-eig
      call power(A,mu,eps,M,eig,conv)
   end if
end if

! If power iteration fails to converge, compute the whole spectrum with Shur decomposition:
if(.not.(conv)) then
  call DGEEV('N','N',n,A,n,wr,wi,vl,1,vr,1,work,lwork,info)
  if(info.eq.0) then
     eig=wr(maxloc(wr,1)) ! Select the largest eigenvalue
     conv=.true.
  else
     print *,'error in DGEEV ',info
  end if
end if

deallocate(work)

return
end subroutine compute_eig

subroutine power(A,mu,eps,M,eig,conv)
use globals
use auxiliary
implicit none
logical :: conv
integer :: i,M,ierr,j
double precision :: A(n,n),v(n),y(n),eig,eps,mu,err,nrm,ndp,mem

conv=.false.

! Make random initial vector
ierr=vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,n,v,0d0,1d0)
v=v/sqrt(sum(v**2))
mem=-100d0

do i=1,M
   call dgemv('N',n,n,1d0,A,n,v,1,0d0,y,1)
   y=y+mu*v
   nrm=sqrt(sum(y**2))
   eig=dot_product(v,y) ! Rayleigh quotient approximates the eigenvalue
   eig=eig-mu
   err=abs(mem-eig)/abs(mem)
   if(err.lt.eps) then
      conv=.true.
      exit
   end if
   v=y/nrm
   mem=eig
end do

return
end subroutine power
