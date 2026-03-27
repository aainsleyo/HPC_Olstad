include 'mkl_vsl.f90'                     ! For quasi-random number generation with the MKL VSL library
module globals
USE MKL_VSL_TYPE                          ! Routines for quasi-random 
USE MKL_VSL                               ! number generation
implicit none
include 'mpif.h'                          ! Fortran MPI header file
integer :: n                              ! Size of the matrices
integer, parameter :: EXIT_TAG=1,DEFAULT_TAG=0,void=0 ! Used to signal workers to return
integer :: proc_num,num_procs             ! Process number and number of processes
end module globals

module auxiliary
use globals
implicit none
TYPE (VSL_STREAM_STATE) :: stream         ! Identifier for the pseudo-random number stream
contains
subroutine get_random_seed(seed)
integer :: un,seed,istat

! Read seed integer from /dev/urandom, a file filled with pseudo-random numbers generated from the state of the computer
open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
   read(un) seed
close(un)

end subroutine get_random_seed

subroutine open_random_stream(seed)
! Setup of stream of pseudo-random numbers using MKL library.
integer :: ierr,seed

! Use a Mersenne Twister pseudorandom number generator; the PRN sequence is identified by the variable "stream"
ierr=vslnewstream(stream,VSL_BRNG_MT19937,seed)
if(ierr.ne.0) then
   print *,'Trouble with the MKL RNG, returns flag ',ierr
end if

end subroutine open_random_stream

subroutine close_random_stream
! Close pseudo-random number stream
integer :: rnd_ierr

rnd_ierr=vsldeletestream( stream )

end subroutine close_random_stream
end module auxiliary

program main
use globals
use auxiliary
implicit none
logical :: ok
double precision :: wtime

call start_MPI(ok)
if(proc_num.eq.0) wtime=MPI_wtime()

if(ok) then
   if(proc_num.eq.0) then
      call manager
      !print *, 'manager called'
   else
      call worker
   end if
end if
if(proc_num.eq.0) then
   wtime=MPI_wtime()-wtime
   print *,'wtime=',wtime,'(s)'
end if
call stop_MPI

end program main

subroutine start_MPI(ok)
use globals
implicit none
! MPI initialization
logical :: my_ok=.true.,ok
integer :: ierr
call mpi_init(ierr)

if(ierr.ne.0) then
   print *,'MPI_init failed!'
   my_ok=.false.
end if

if(my_ok) then
   call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)
   if(ierr.ne.0) then
      print *,'MPI_comm_size failed!'
      my_ok=.false.
   end if
   if(my_ok) then
      call mpi_comm_rank(MPI_COMM_WORLD,proc_num,ierr)
      if(ierr.ne.0) then
         print *,'MPI_comm_rank failed!'
         my_ok=.false.
      end if
   end if
end if
! Check if everyone is ok.
call mpi_allreduce(my_ok,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
end subroutine start_MPI
subroutine stop_MPI
use globals
implicit none
logical :: init
integer :: ierr

! Wait until everybody is done. One process "finalize"ing while others are still working can cause ugly crashes.
call mpi_barrier(MPI_COMM_WORLD,ierr)
! Check if MPI has been initialized.
call mpi_initialized(init,ierr)
! If it is, call finalize.
if(init) call mpi_finalize(ierr)

end subroutine stop_MPI
