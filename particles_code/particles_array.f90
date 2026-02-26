! Global variables
module globals
implicit none
integer :: n=100                               ! number of particles
double precision :: L=10.0d0
double precision :: rc
integer :: M                                   ! number of sectors per side
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant

contains

subroutine update_M()
    M = floor(L/rc)
end subroutine update_M

end module globals

!module griding
module griding
implicit none
contains

function TWOtoONE(x, y, box_size) result(id)
    integer, intent(in) :: x, y, box_size
    integer :: id
    id = y * box_size + x
end function TWOtoONE

end module griding

!module sectors
module sectors
use globals
use griding
implicit none
contains

integer function indexfxn(px, py) result(idx)

    double precision, intent(in) :: px, py
    integer :: ix, iy

    ix = floor((px + L/2d0) / rc)
    iy = floor((py + L/2d0) / rc)

    if (ix < 0 .or. ix >= M .or. iy < 0 .or. iy >= M) then
        print *, "Sector error"
        stop
    end if

    idx = TWOtoONE(ix, iy, M)

end function indexfxn

end module sectors

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
logical, allocatable, dimension(:) :: is_tracked
double precision :: dt,kT,g,mass,sigma,eps      ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0 ! particle positions, accellerations, velocities, half-step velocities, initial positions
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.0001d0 ! time step size
kT=1d0    ! energy
g=1d0     ! drag coefficient
mass=1d0     ! mass of the particles, can be normalized to 1.
sigma=1d-3              ! Potential parameters
eps=1d0
rc= 2.0d0 !sigma*2d0**(1d0/6d0) ! Effective particle size
call update_M()

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
   gr1=sqrt(kT/(mass))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(mass))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
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

!module reys here?
!module sectors here?

program main
use globals
use Langevin
use BC
use sectors
implicit none
integer :: i,j
double precision :: t,t_max,m1,m2,rx,ry,dij,F
double precision :: wtime,begin,end
double precision, allocatable, dimension(:) :: ran1,ran2

! Open files
open(11,file='trajectories')
open(12,file='means')
open(15,file='velocities')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),is_tracked(n),ran1(n),ran2(n))

is_tracked = .True.
t=0d0
t_max=0.1d0     ! integration time

call set_parameters
call initialize_particles

print *, "L =", L
print *, "rc =", rc
print *, "M =", M

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

   do i=1,n
      if(abs(x(i)) > L/2 .or. abs(y(i)) > L/2) then
        print *, "BOUNDARY ERROR"
        stop
      end if
   end do

      
   ax=0d0                   ! Add forces here if any
   ay=0d0                   ! Add forces here if any

   call random_number(ran1)
   ran1=ran1-0.5d0
   call random_number(ran2)
   ran2=ran2-0.5d0
      
   ax=ax-pref1*vhx+pref2*ran1
   ay=ay-pref1*vhy+pref2*ran2

   !$omp parallel do private(j,rx,ry,dij,F)
   do i=1,n-1
      do j=i+1,n
         if(j.ne.i) then
            rx=x(j)-x(i)
            ry=y(j)-y(i)
            dij=sqrt(rx**2 + ry**2)
            if(dij.lt.rc) then
               !print *,'interaction detected at t=',t,' (x,y)=',x(i),y(i)
               F=4d0*eps*( -12d0*sigma**12/dij**13 + 6D0* sigma**6/dij**7 )
               
               !$omp atomic
               ax(i)=ax(i)+F*rx/(dij*mass)
               
               !$omp atomic
               ay(i)=ay(i)+F*ry/(dij*mass)
            end if
         end if
      end do
   end do
   !$omp end parallel do
   vx=vhx+0.5d0*ax*dt
   vy=vhy+0.5d0*ay*dt

   do i=1,n
      write(11,*) x(i),y(i)
   enddo

   do i=1,n
      write(15,*) vx(i),vy(i)
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
close(15)

end program main
