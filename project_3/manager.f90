subroutine manager
use globals
use auxiliary
implicit none
logical :: start(1:num_procs-1)
integer :: seed,ierr,recvd=0,worker,exit_tags_sent,tag,ndat,p,failed
integer, dimension(mpi_status_size) :: status
double precision :: buf
double precision, allocatable, dimension(:) :: eigs

! Set matrix size and number of eigenvalues
call init(ndat)
! Broadcast matrix size
call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Allocate array for received eigenvalues
allocate(eigs(ndat))

start=.true.
exit_tags_sent=0
failed=0
open(22,file='eigs')

! Receive data until at least ndat have been received
do while(.true.)
   ! For each worker check if a seed is needed
   do p=1,num_procs-1
      if(start(p)) then
         call get_random_seed(seed)
         if(recvd.lt.ndat-num_procs+2) then
            tag=DEFAULT_TAG
         else
            tag=EXIT_TAG
            exit_tags_sent=exit_tags_sent+1
         end if
         call mpi_send(seed,1,MPI_INTEGER,p,tag,MPI_COMM_WORLD,ierr)
!         print *,'Manager sent seed ',seed,' to worker ',p
         start(p)=.false.
      end if
   end do
   if(exit_tags_sent.eq.num_procs-1) then
      print *,'Manager exiting...'
      exit
   end if
   
   call mpi_recv(buf,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
   worker = status(MPI_SOURCE)
   tag = status(MPI_TAG)
   if(tag.eq.0) then
      recvd=recvd+1
      eigs(recvd)=buf
      !print *. "eig"
      write(22,'(i4,e14.7)') worker,buf
   else
      failed=failed+1
      print '(i5,a7,i5,a8)',recvd,' received, ',failed,' failed'
   end if
   start(worker)=.true.
end do

close(22)

deallocate(eigs)
return
end subroutine manager

subroutine init(ndat)
use globals
implicit none
integer :: ndat

! Matrix dimension and sample size.. hard-coded now, but can be read from command line or file...
!read(*,*) n
!read(*,*) ndat

n=500
ndat=1000

return
end subroutine init
