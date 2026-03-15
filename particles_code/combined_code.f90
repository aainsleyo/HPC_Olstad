module globals
! Global variables
use omp_lib                                    ! help the compiler find the OMP libraries
implicit none
integer :: n=100                               ! number of particles
double precision :: L=1.0d0
double precision, parameter :: pi=2.0d0*asin(1.0d0) ! numerical constant
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
kT=1d0    ! energy
g=1d0     ! drag coefficient
m=1d0     ! mass of the particles, can be normalized to 1.
sigma= 0.01d0 !1d-3              ! Potential parameters
eps=1d0
rc=sigma*2d0**(1d0/6d0) ! Effective particle size

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i, ix, iy, idx
integer :: ngrid
double precision :: spacing
double precision :: ran1(n),ran2(n),gr1(n),gr2(n)
! Give particles initial position and velocity

!Grid Init
   ngrid   = ceiling(sqrt(real(n, 8)))           ! e.g. n=1000 -> ngrid=32
   spacing = L / real(ngrid, 8)                  ! gap between grid points

   idx = 1
   outer: do iy = 0, ngrid-1
     do ix = 0, ngrid-1
        if (idx .gt. n) exit outer
        x(idx) = -L/2.0d0 + (real(ix, 8) + 0.5d0) * spacing
        y(idx) = -L/2.0d0 + (real(iy, 8) + 0.5d0) * spacing
        idx = idx + 1
     end do
    end do outer

   call random_seed()
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
  integer, parameter :: b=10
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
     do while (y(i) .GT. L/2)
         y(i) = L - y(i) 
         vhy(i) = -vhy(i)
     end do

     !> Left boundary
     do while(x(i) .LT. -L/2)
         x(i) = -L - x(i) 
         vhx(i) = -vhx(i)
     end do

     !> Bottom boundary
     do while(y(i) .LT. -L/2)
         y(i) = -L - y(i)
         vhy(i) = -vhy(i)
     end do

     !> Right boundary
     do while(x(i) .GT. L/2)
         x(i) = L - x(i)
         vhx(i) = -vhx(i)
     end do

     !> Final check
     if(abs(x(i)).GT.L/2 .OR. abs(y(i)).GT.L/2) then
         !> Particle is still outside, don't track it
         is_tracked(i) = .FALSE.
         !write(*,*) 'WARNING', x(i), y(i)
     end if
end subroutine impose_BC
end module BC

program main
use globals
use domainDecomposition
use Langevin
use BC
implicit none
integer :: i,j,lim(0:b*b,2),s,ns,p1,p2,step
double precision :: t,t_max,m1,m2,rx,ry,dij,F
double precision :: wtime,t_begin,t_end,Tint,TsinceWrite
double precision, allocatable, dimension(:) :: ran1,ran2, scrapx, scrapy, vxscrap, vyscrap

! Open files
open(11,file='trajectories')
open(12,file='means')
open(13,file='testing')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),is_tracked(n),ran1(n),ran2(n), scrapx(n), scrapy(n), vxscrap(n), vyscrap(n))

call buildNBL()

is_tracked = .True.
t=0d0
t_max=1.0d0     ! integration time
Tint = 0.01d0
TsinceWrite = 0d0
step = 0

call set_parameters
call initialize_particles

t_begin = omp_get_wtime()
!call cpu_time(begin)

! Conclusion: we need to re-order the loop like this:
! a. update half-step velocities
! b. update positions
! c. compute accellerations/forces
! d. update all velocities

!$omp parallel private(i,j,s,ns,p1,p2,rx,ry,dij,F)
do while(t.lt.t_max)
   ! one thread: writing to disk
   ! one thread: fetch psuedo-random numbers
   ! one thread: update velocity, position, impose Bc

   !$omp single
   scrapx=x
   scrapy=y
   vxscrap=vx
   vyscrap=vy
   !$omp end single

   !$omp sections
   !$omp section
   ! Thread 1
   if(tSinceWrite.gt.Tint) then
   do i=1,n
      write(11,*) scrapx(i), scrapy(i)
   end do
   write(12,*) t, sum(m*(vxscrap**2+vyscrap**2)/(2*n))
   TsinceWrite = 0d0
   end if

   !$omp section
   ! Thread 2
   call random_number(ran1)
   ran1 = ran1 - 0.5d0
   call random_number(ran2)
   ran2 = ran2 - 0.5d0

   !$omp section
   ! Thread 3
   vhx = vx + 0.5d0*ax*dt
   vhy = vy + 0.5d0*ay*dt
   x   = x  + vhx*dt
   y   = y  + vhy*dt

   do j=1,n
      call impose_BC(j)
   end do

   call order(x,y,vx,vy,x0,y0,vhx,vhy,lim)
   
   ax=0d0      !add any forces here
   ay=0d0      !add any forces here

   !$omp end sections
   !$omp flush
   !$omp single

   ax = -pref1*vhx + pref2*ran1
   ay = -pref1*vhy + pref2*ran2
   !$omp end single
   !$omp do private(ns,p1,p2,rx,ry,dij,F)
   do s=0,b*b-1
      do i=1,9
         if(nbl(s,ns).eq.-1) exit
         ns=nbl(s,i)
         do p1=lim(s,1),lim(s,2)
            do p2=lim(nbl(s,ns),1),lim(nbl(s,ns),2)
               if(p1.eq.p2) cycle
               rx=x(p2)-x(p1)
               ry=y(p2)-y(p1)
               dij=sqrt(rx**2+ry**2)
               if(dij.lt.rc) then
                  F=4d0*eps*(-12d0*sigma**12/dij**13 &
                            + 6d0*sigma**6/dij**7)
                  !$omp atomic
                  ax(p1)=ax(p1)+F*rx/(dij*m)
                  !$omp atomic
                  ay(p1)=ay(p1)+F*ry/(dij*m)
               end if
            end do
         end do
      end do
   end do
   !$omp end do

   !$omp single
   vx = vhx + 0.5d0*ax*dt
   vy = vhy + 0.5d0*ay*dt
   t = t + dt 
   step = step + 1
   tSinceWrite=tSinceWrite + dt
   !$omp end single

end do
!$omp end parallel

t_end = omp_get_wtime()
!call cpu_time(end)

print *,'Wtime=',t_end - t_begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,x0,y0,is_tracked,ran1,ran2, scrapx, scrapy, vxscrap, vyscrap)
! Close files
close(11)
close(12)
close(13)

end program main
