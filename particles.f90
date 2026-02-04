module globals
! Global variables
implicit none
integer :: n=1000                               ! number of particles
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
double precision, parameter :: L=1.0d0
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
logical, allocatable, dimension(:) :: outside
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0 ! particle positions, accellerations, velocities, half-step velocities, initial positions
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.01d0 ! time step size
kT=1d0    ! energy
g=1d0     ! drag coefficient
m=1d0     ! mass of the particles, can be normalized to 1.

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2
! Give particles initial position and velocity
do i=1,n
   call random_number(ran1)                       ! uses the built-in PRNG, easy but not very accurate
   call random_number(ran2)
   x(i)=L*(ran1-0.5d0)
   x0(i)=x(i)
   y(i)=L*(ran2-0.5d0)
   y0(i)=y(i)
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
end do

end subroutine initialize_particles
end module Langevin

module BC
  ! Subroutines related to the boundary conditions
  use globals
  use Langevin
  implicit none
contains
  subroutine impose_BC(i)
    integer :: i

    ! decide on BC and implement those here

    check if (x(i), y(i)) is inside [0,L] X [0,L]
    if NOT:
       check top: if y(i) > L/2 then:
          reflect y(i) in y=L/2
          reflect vhy(i)
       check left: if x(i) < -L/2 then:
	  reflect x(i) in x=-L/2
	  reflect vhx(i)
       check bottom: if y(i) < -L/2 then:
          reflect y(i) in y=-L/2
          reflect vhy(i)
       check left: if x(i) > L/2 then:
          reflect x(i) in x=L/2
          reflect vhx(i)
    if particle i is still outside then set outside(i) = True


  end subroutine impose_BC
end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i
double precision :: t,t_max,ran1,ran2,m1,m2
double precision :: wtime,begin,end

! Open files
open(12,file='means')

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),outside(n))

outside = .False.
t=0d0
t_max=100d0     ! integration time

call set_parameters
call initialize_particles

call cpu_time(begin)

do while(t.lt.t_max)
   do i=1,n     ! run the velocity-Verlet algorithm...
      if(outside(i)) continue
      vhx(i)=vx(i)+0.5d0*ax(i)*dt
      vhy(i)=vy(i)+0.5d0*ay(i)*dt
      x(i)=x(i)+vhx(i)*dt
      y(i)=y(i)+vhy(i)*dt

      do j=1,n
	call impose_BC(j)
      end do
      
      ax(i)=0d0                   ! Add forces here if any
      ay(i)=0d0                   ! Add forces here if any

      call random_number(ran1)
      ran1=ran1-0.5d0
      call random_number(ran2)
      ran2=ran2-0.5d0
      ax(i)=ax(i)-pref1*vhx(i)+pref2*ran1
      ay(i)=ay(i)-pref1*vhy(i)+pref2*ran2
      
      vx(i)=vhx(i)+0.5d0*ax(i)*dt
      vy(i)=vhy(i)+0.5d0*ay(i)*dt
   end do
   t=t+dt
   write(12,*) t,sqrt(sum((x-x0)**2+(y-y0)**2)/real(n,8))
end do

call cpu_time(end)
print *,'Wtime=',end-begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,x0,y0,outside)
! Close files
close(11)
close(12)

end program main
