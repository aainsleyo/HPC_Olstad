module globals
! Global variables
use omp_lib                                    ! help the compiler find the OMP libraries
implicit none
integer :: n=1000                               ! number of particles
double precision, parameter :: L=1.0d0
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
logical, allocatable, dimension(:) :: is_tracked
double precision :: dt,kT,g,m,sigma,eps,rc      ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0 ! particle positions, accellerations, velocities, half-step velocities, initial positions
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.00001d0 ! time step size
kT=0.1d0    ! energy
g=0.1d0     ! drag coefficient
m=1d0     ! mass of the particles, can be normalized to 1.
sigma=1d-3              ! Potential parameters
eps=1d0
rc=sigma*2d0**(1d0/6d0) ! Effective particle size

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i
double precision :: ran1(n),ran2(n),gr1(n),gr2(n)
! Give particles initial position and velocity

   call random_seed(put=(/1/))
   call random_number(ran1)                       ! uses the built-in PRNG, easy but not very accurate
   call random_number(ran2)
   x=L*(ran1-0.5d0)
   x0=x
   y=L*(ran2-0.5d0)
   y0=y
   ax=0d0
   ay=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx=gr1
   vy=gr2

end subroutine initialize_particles
end module Langevin

module domainDecomposition
  use globals
  use Langevin
  implicit none
  integer, parameter :: b=4
  integer :: nbl(0:b*b-1,9)
contains
  include "neighbourlist.f90"
  include "putParticleInBox.f90"
  include "sortParticles.f90"
end module domainDecomposition

module BC
   ! Subroutines related to the boundary conditions
   use globals
   use langevin
   implicit none
contains
   subroutine impose_BC(i)
     integer :: i
     !> Recall we are inside an LxL box centered at the origin

     !> Top boundary
     if(y(i) .GT. L/2) then
         y(i) = L - y(i) 
         vhy(i) = -vhy(i)
     end if

     !> Left boundary
     if(x(i) .LT. -L/2) then
         x(i) = -L - x(i) 
         vhx(i) = -vhx(i)
     end if

     !> Bottom boundary
     if(y(i) .LT. -L/2) then
         y(i) = -L - y(i)
         vhy(i) = -vhy(i)
     end if

     !> Right boundary
     if(x(i) .GT. L/2) then
         x(i) = L - x(i)
         vhx(i) = -vhx(i)
     end if

   end subroutine impose_BC
end module BC

program main
  use globals
  use domainDecomposition
use Langevin
use BC
implicit none
integer :: i,j,lim(0:b*b,2),s,ns,p1,p2
double precision :: t,t_max,m1,m2,rx,ry,dij,F
double precision :: wtime,begin,end,Tint,TsinceWrite
double precision, allocatable, dimension(:) :: ran1,ran2,xscrap,yscrap,vxscrap,vyscrap

! Open files
open(11,file='trajectories')
open(12,file='means')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),is_tracked(n),ran1(n),ran2(n),xscrap(n),yscrap(n),vxscrap(n),vyscrap(n))

call buildNBL()

is_tracked = .True.
t=0d0
t_max=1.0d0     ! integration time
Tint = 0.001d0
TsinceWrite=0d0
call set_parameters
call initialize_particles

begin = omp_get_wtime()
!call cpu_time(begin)

! Conclusion: we need to re-order the loop like this:
! a. update half-step velocities
! b. update positions
! c. compute accellerations/forces
! d. update all velocities

!$omp parallel
do while(t.lt.t_max)
   ! one thread: fetch pseudo-random numbers
   ! one thread: update velocity, position, impose BC
   xscrap=x
   yscrap=y
   vxscrap=vx
   vyscrap=vy
   !$omp sections
   !$omp section  ! Write to disk
      if(tSinceWrite.gt.Tint) then
         do i=1,n
            write(11,*) xscrap(i),yscrap(i)
         enddo
!         write(11,*) ' '
         write(12,*) t,sqrt(sum((xscrap-x0)**2+(yscrap-y0)**2)/real(n,8))
         tSinceWrite=0d0
      end if
   !  write(12,*) t, sum(m*(vxscrap**2+vyscrap**2)/(2*n))
   !$omp section
      call random_number(ran1)
      ran1=ran1-0.5d0
      call random_number(ran2)
      ran2=ran2-0.5d0
   !$omp section
      vhx=vx+0.5d0*ax*dt
      vhy=vy+0.5d0*ay*dt
      x=x+vhx*dt
      y=y+vhy*dt
      do j=1,n
         call impose_BC(j)
      end do
      call order(x,y,vx,vy,x0,y0,vhx,vhy,lim) ! BUGFIX 11/03: half-step velocities must also be re-ordered.
      ax=0d0                   ! Add forces here if any
      ay=0d0                   ! Add forces here if any
   !$omp end sections
   !$omp flush
   !$omp single
      ax=ax-pref1*vhx+pref2*ran1
      ay=ay-pref1*vhy+pref2*ran2
   !$omp end single

! Our first attempt at parallelization of the code: run the computation of the distances and interaction forces on multiple threads:
!$omp do private(s,i,ns,p1,p2,rx,ry,dij,F)
   do s=0,b*b-1
      do i=1,9
         if(nbl(s,i).eq.-1) exit
         ns=nbl(s,i)  ! BUGFIX 11/03: the loop counter was used as sector index instead of the entries of nbl.
         do p1=lim(s,1),lim(s,2)
            do p2=lim(ns,1),lim(ns,2)
               if(p1.eq.p2) cycle
               rx=x(p2)-x(p1)
               ry=y(p2)-y(p1)
               dij=sqrt(rx**2 + ry**2)
               if(dij.lt.rc) then
                  F=4d0*eps*( -12d0*sigma**12/dij**13 + 6D0* sigma**6/dij**7 )
                  ax(p1)=ax(p1)+F*rx/(dij*m)
                  ay(p1)=ay(p1)+F*ry/(dij*m)
               end if
            end do
         end do
      end do
   end do
   !$omp end do
   !$omp single
      vx=vhx+0.5d0*ax*dt
      vy=vhy+0.5d0*ay*dt
      t=t+dt
      tSinceWrite=tSinceWrite+dt
   !$omp end single
end do
!$omp end parallel

end = omp_get_wtime()
!call cpu_time(end)
print *,'Wtime=',end-begin
if(lim(b*b,2).gt.lim(b*b,1)) print '(a5,i7,a8,i8,a11)','Lost ',lim(b*b,2)-lim(b*b,1),' out of ',n,' particles.'
! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,x0,y0,is_tracked,ran1,ran2,xscrap,yscrap,vxscrap,vyscrap)
! Close files
close(11)
close(12)

end program main
