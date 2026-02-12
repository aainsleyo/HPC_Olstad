module globals
! Global variables
implicit none
integer :: n=100                               ! number of particles
double precision :: L=1.0d0
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
dt=0.0001d0 ! time step size
kT=1d0    ! energy
g=1d0     ! drag coefficient
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

     !> Final check
     if(abs(x(i)).GT.L/2 .OR. abs(y(i)).GT.L/2) then
         !> Particle is still outside, don't track it
         is_tracked(i) = .FALSE.
     endif

   end subroutine impose_BC
end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i,j
double precision :: t,t_max,m1,m2,rx,ry,dij,F
double precision :: wtime,begin,end
double precision, allocatable, dimension(:) :: ran1,ran2

! Open files
open(11,file='trajectories')
open(12,file='means')
open(13,file='testing')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),is_tracked(n),ran1(n),ran2(n))

is_tracked = .True.
t=0d0
t_max=0.1d0     ! integration time

call set_parameters
call initialize_particles

call cpu_time(begin)

! Conclusion: we need to re-order the loop like this:
! a. update half-step velocities
! b. update positions
! c. compute accellerations/forces
! d. update all velocities
do while(t.lt.t_max)
   vhx=vx+0.5d0*ax*dt
   vhy=vy+0.5d0*ay*dt
   x=x+vhx*dt
   y=y+vhy*dt
   do j=1,n
      call impose_BC(j)
   end do
      
   ax=0d0                   ! Add forces here if any
   ay=0d0                   ! Add forces here if any

   call random_number(ran1)
   ran1=ran1-0.5d0
   call random_number(ran2)
   ran2=ran2-0.5d0
      
   ax=ax-pref1*vhx+pref2*ran1
   ay=ay-pref1*vhy+pref2*ran2

   do i=1,n
      do j=1,n
         if(j.ne.i) then
            rx=x(j)-x(i)
            ry=y(j)-y(i)
            dij=sqrt(rx**2 + ry**2)
            if(dij.lt.rc) then
               write(13,*) x(i),y(i)
               print *,'interaction detected at t=',t,' (x,y)=',x(i),y(i)
               F=4d0*eps*( -12d0*sigma**12/dij**13 + 6D0* sigma**6/dij**7 )
               ax(i)=ax(i)+F*rx/(dij*m)
               ay(i)=ay(i)+F*ry/(dij*m)
            end if
         end if
      end do
   end do
   vx=vhx+0.5d0*ax*dt
   vy=vhy+0.5d0*ay*dt

   do i=1,n
      write(11,*) x(i),y(i)
   enddo

   t=t+dt
   write(12,*) t,sqrt(sum((x-x0)**2+(y-y0)**2)/real(n,8))
end do

call cpu_time(end)
print *,'Wtime=',end-begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,x0,y0,is_tracked,ran1,ran2)
! Close files
close(11)
close(12)
close(13)

end program main
