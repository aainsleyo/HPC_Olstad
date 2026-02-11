module globals
! Global variables
implicit none
integer :: n=100                               ! number of particles
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
double precision, parameter :: L=1.0d0
end module globals

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
logical, allocatable, dimension(:) :: outside
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: sigma,eps,rc                ! additional parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision :: x_trial, y_trial
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0 ! particle positions, accellerations, velocities, half-step velocities, initial positions
contains

subroutine set_parameters

! Set time step and physical parameters
dt=0.01d0 ! time step size
kT=1d0    ! energy
g=1d0     ! drag coefficient
m=1d0     ! mass of the particles, can be normalized to 1.

sigma=0.05d0     ! potentil parameters (sigma, eps, rc)
eps=1d0
rc=sigma*24d8**(1d0/6d0) ! effective particle size

! Set auxiliary parameters
pref1=g
pref2=sqrt(2d0*kT*g/dt)
!pref2=sqrt(24d0*kT*g/dt)
end subroutine set_parameters

subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2

! Give particles initial position and velocity
do i = 1,n
	call random_number(ran1)                       ! uses the built-in PRNG, easy but not very accurate
	call random_number(ran2)

	x(i)=L*(ran1-0.5d0)
	y(i)=L*(ran2-0.5d0)
	
	x0(i) = x(i)
	y0(i) = y(i)

	ax=0d0
	ay=0d0

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
    integer, intent(in) :: i
    logical :: still_outside

    still_outside = .false.

    ! --- Top wall (y = +L/2)
    if (y(i) > 0.5d0*L) then
       y(i)  =  L - y(i)           ! reflect position
       vhy(i) = -vhy(i)            ! reflect velocity
    end if

    ! --- Bottom wall (y = -L/2)
    if (y(i) < -0.5d0*L) then
       y(i)  = -L - y(i)
       vhy(i) = -vhy(i)
    end if

    ! --- Right wall (x = +L/2)
    if (x(i) > 0.5d0*L) then
       x(i)  =  L - x(i)
       vhx(i) = -vhx(i)
    end if

    ! --- Left wall (x = -L/2)
    if (x(i) < -0.5d0*L) then
       x(i)  = -L - x(i)
       vhx(i) = -vhx(i)
    end if

    ! --- Final safety check
    if ( x(i) >  0.5d0*L .or. x(i) < -0.5d0*L .or. &
         y(i) >  0.5d0*L .or. y(i) < -0.5d0*L ) then
       outside(i) = .true.
    end if

  end subroutine impose_BC

end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i, j
double precision :: t,t_max,ran1,ran2,m1,m2, dij, rx, ry, F
double precision :: wtime,begin,end

! Open files
open(11,file="positions.dat",status="replace")
open(12,file='means')
open(14,file="interactions.dat",status="replace")

! Open trajectory files (for debugging / testing)
open(20, file="traj_1.dat", status="replace")
open(21, file="traj_2.dat", status="replace")
open(22, file="traj_3.dat", status="replace")

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),outside(n))

outside = .False.
t=0d0
t_max=5.0d0     ! integration time

call set_parameters
call initialize_particles

call cpu_time(begin)

! conclusion we need ot re order the loop like this:
! a. update half step velocities
! b. update positions 
! c. compute accellerations
! d. update all velocities
! now we can compare the positions to whats happening now and not the time step in the past...? I think.
do while (t < t_max)
       
   ! --- Half-step velocities
   vhx = vx + 0.5d0*ax*dt
   vhy = vy + 0.5d0*ay*dt

   ! --- Position Update
   x = x + vhx*dt
   y = y + vhy*dt

      ! ======================================
      ! TEST 2: OVERSHOOT TEST (BEFORE APPLYING BC)
      ! ======================================
      !x_trial = x(i) + vhx(i)*dt
      !y_trial = y(i) + vhy(i)*dt

      !if (abs(x_trial) > 0.5d0*L .or. abs(y_trial) > 0.5d0*L) then
      !   print *, "OVERSHOOT DETECTED:"
      !   print *, " i =", i
      !   print *, " t =", t
      !   print *, " x_trial,y_trial =", x_trial, y_trial
      !   print *, " vhx,vhy =", vhx(i), vhy(i)
      !end if
      ! ======================================

   ! --- impose boundary conditions
   do i = 1,n
      call impose_BC(i)
   end do

   ! =====================================================
   ! TEST 1: particle must be inside box or flagged outside
   ! =====================================================
   do i = 1,n
      if (.not. outside(i)) then
         if (abs(x(i)) > 0.5d0*L .or. abs(y(i)) > 0.5d0*L) then
            print *, "BC ERROR:"
            print *, " i =", i
            print *, " t =", t
            print *, " x,y =", x(i), y(i)
            print *, " vhx,vhy =", vhx(i), vhy(i)
            stop
         end if
      end if
   end do
   ! =====================================================
      
      ! --- reset accelerations
      ax = 0d0
      ay = 0d0

      ! --- Pair interactions
      do i = 1,n
         do j = i,n
             rx=x(j)-x(i)
             ry=y(j)-y(i)
             dij=sqrt(rx*rx + ry*ry)

             if(dij.lt.rc) then
               F=4d0*eps*(-12d0*sigma**12/dij**13 + 6D0*sigma**6/dij**7)
               write(14,*) x(i), y(i)
               print *, "interaction detected at t=",x(i), y(i)
               ax(i)=ax(i)+F*rx/(dij*m)
               ay(i)=ay(i)+F*ry/(dij*m)

               !ax(j)=ax(j)-F*rx/(dij*m)
               !ay(j)=ay(j)-F*ry/(dij*m)

             end if
        end do
      end do

   ! --- Langevin forces (drag + noise)
   do i=1,n
      
      call random_number(ran1)
      ran1 = ran1 - 0.5d0
      call random_number(ran2)
      ran2 = ran2 - 0.5d0

      ax(i) = ax(i) - pref1*vhx(i) + pref2*ran1
      ay(i) = ay(i) - pref1*vhy(i) + pref2*ran2


   end do

   ! --- full-step velocities
   vx = vhx + 0.5d0*ax*dt
   vy = vhy + 0.5d0*ay*dt

   ! --- Diagnostics
   write(11,*) (x(i),y(i), i=1,n)
   write(12,*) t, sqrt(sum((x-x0)**2+(y-y0)**2)/real(n,8))

   if (n>=1) write(20,*) t,x(1),y(1)
   if (n>=2) write(21,*) t,x(2),y(2)
   if (n>=3) write(22,*) t,x(3),y(3)  

   t = t + dt
end do

call cpu_time(end)
print *,'Wtime=',end-begin

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,vhx,vhy,x0,y0,outside)
! Close files
close(11)
close(12)
close(14)
close(20)
close(21)
close(22)

end program main
