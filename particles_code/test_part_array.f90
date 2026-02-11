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
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0
contains

subroutine set_parameters

! Set time step and physical parameters
dt=0.01d0 ! time step size
kT=1d0    ! energy
g=1d0     ! drag coefficient
m=1d0     ! mass of the particles, can be normalized to 1.

! CORRECTED: Use reasonable sigma for interactions
sigma=1d-3    
eps=1d0
rc=sigma*24d0**(1d0/6d0)

! Set auxiliary parameters
pref1=g
pref2=sqrt(2d0*kT*g/dt) 

print *, "Parameters set:"
print *, "  sigma =", sigma
print *, "  rc =", rc
print *, "  pref2 =", pref2
print *, ""

end subroutine set_parameters

subroutine initialize_particles
integer :: i
double precision :: ran1,ran2,gr1,gr2

! Give particles initial position and velocity
do i = 1,n
	call random_number(ran1)
	call random_number(ran2)

	x(i)=L*(ran1-0.5d0)
	y(i)=L*(ran2-0.5d0)
	
	x0(i) = x(i)
	y0(i) = y(i)

	call random_number(ran1)
	call random_number(ran2)

	gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
	gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
	
	vx(i)=gr1
	vy(i)=gr2

end do

! Initialize accelerations
ax=0d0
ay=0d0

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

    ! --- Safety check: try one more reflection if still outside
    if ( x(i) >  0.5d0*L .or. x(i) < -0.5d0*L .or. &
         y(i) >  0.5d0*L .or. y(i) < -0.5d0*L ) then
       
       if (y(i) > 0.5d0*L) then
          y(i) = L - y(i)
          vhy(i) = -vhy(i)
       end if
       if (y(i) < -0.5d0*L) then
          y(i) = -L - y(i)
          vhy(i) = -vhy(i)
       end if
       if (x(i) > 0.5d0*L) then
          x(i) = L - x(i)
          vhx(i) = -vhx(i)
       end if
       if (x(i) < -0.5d0*L) then
          x(i) = -L - x(i)
          vhx(i) = -vhx(i)
       end if
       
       if ( x(i) >  0.5d0*L .or. x(i) < -0.5d0*L .or. &
            y(i) >  0.5d0*L .or. y(i) < -0.5d0*L ) then
          outside(i) = .true.
          print *, "WARNING: Particle",i,"escaped at x=",x(i),"y=",y(i)
       end if
    end if

  end subroutine impose_BC

end module BC

program main
use globals
use Langevin
use BC
implicit none
integer :: i, j, interaction_count, total_interactions
double precision :: t,t_max,ran1,ran2,gr1,gr2, dij, rx, ry, F
double precision :: begin,end

! Open files
open(11,file="positions.dat",status="replace")
open(12,file='means.dat',status="replace")
open(13,file='interactions.dat',status="replace")

! Open trajectory files (for debugging / testing)
open(20, file="traj_1.dat", status="replace")
open(21, file="traj_2.dat", status="replace")
open(22, file="traj_3.dat", status="replace")

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n),outside(n))

outside = .False.
t=0d0
t_max=5.0d0     ! integration time
total_interactions = 0

call set_parameters
call initialize_particles

print *, "Starting simulation:"
print *, "  n =", n, "particles"
print *, "  t_max =", t_max
print *, "  dt =", dt
print *, ""

call cpu_time(begin)

! VELOCITY VERLET INTEGRATION (array-based, physically correct)
do while (t < t_max)
       
   ! ===================================================
   ! STEP 1: Half-step velocities (using old accelerations)
   ! ===================================================
   vhx = vx + 0.5d0*ax*dt
   vhy = vy + 0.5d0*ay*dt

   ! ===================================================
   ! STEP 2: Position update (all particles together)
   ! ===================================================
   x = x + vhx*dt
   y = y + vhy*dt

   ! ===================================================
   ! STEP 3: Apply boundary conditions to all particles
   ! ===================================================
   do i = 1,n
      call impose_BC(i)
   end do

   ! ===================================================
   ! STEP 4: Sanity check
   ! ===================================================
   do i = 1,n
      if (.not. outside(i)) then
         if (abs(x(i)) > 0.5d0*L .or. abs(y(i)) > 0.5d0*L) then
            print *, "BC ERROR: particle",i,"at x=",x(i),"y=",y(i)
            stop
         end if
      end if
   end do

   ! ===================================================
   ! STEP 5: Reset accelerations
   ! ===================================================
   ax = 0d0
   ay = 0d0

   ! ===================================================
   ! STEP 6: Compute pair interactions (CORRECT - no double counting)
   ! ===================================================
   interaction_count = 0
   do i = 1,n-1
      do j = i+1,n  ! Each pair computed ONCE

          rx=x(j)-x(i)
          ry=y(j)-y(i)
          dij=sqrt(rx*rx+ry*ry)

          if(dij.lt.rc) then
             interaction_count = interaction_count + 1
             total_interactions = total_interactions + 1

             ! Lennard-Jones force
             F=4d0*eps*(-12d0*sigma**12/dij**13+6d0*sigma**6/dij**7)

             ! Apply force to particle i
             ax(i)=ax(i)+F*rx/(dij*m)
             ay(i)=ay(i)+F*ry/(dij*m)

             ! Apply equal and opposite force to particle j (Newton's 3rd law)
             ax(j)=ax(j)-F*rx/(dij*m)
             ay(j)=ay(j)-F*ry/(dij*m)

             ! Log interaction
             write(13,'(A,F8.4,A,I4,A,I4,A,ES12.4,A,ES12.4)') &
                  't=', t, ' particles ', i, '-', j, &
                  ' dij=', dij, ' F=', F

          end if
     end do
   end do

   ! Print interaction summary
   if (interaction_count > 0) then
      print '(A,F8.4,A,I5,A)', 't=', t, ': ', interaction_count, ' interactions'
   end if

   ! ===================================================
   ! STEP 7: Add Langevin forces (drag + Gaussian noise)
   ! ===================================================
   do i=1,n
      
      ! Generate Gaussian random numbers (Box-Mueller)
      call random_number(ran1)
      call random_number(ran2)
      gr1 = sqrt(-2d0*log(ran1))*cos(2d0*pi*ran2)
      gr2 = sqrt(-2d0*log(ran1))*sin(2d0*pi*ran2)

      ! Drag force and thermal noise
      ax(i) = ax(i) - pref1*vhx(i) + pref2*gr1
      ay(i) = ay(i) - pref1*vhy(i) + pref2*gr2

   end do

   ! ===================================================
   ! STEP 8: Full-step velocities (using new accelerations)
   ! ===================================================
   vx = vhx + 0.5d0*ax*dt
   vy = vhy + 0.5d0*ay*dt

   ! ===================================================
   ! STEP 9: Output diagnostics
   ! ===================================================
   write(11,*) (x(i),y(i), i=1,n)
   write(12,*) t, sqrt(sum((x-x0)**2+(y-y0)**2)/real(n,8))

   if (n>=1) write(20,*) t,x(1),y(1)
   if (n>=2) write(21,*) t,x(2),y(2)
   if (n>=3) write(22,*) t,x(3),y(3)  

   t = t + dt
end do

call cpu_time(end)

print *, ""
print *, "======================"
print *, "Simulation complete!"
print *, "======================"
print *, "Wall time: ", end-begin, " seconds"
print *, "Total interactions: ", total_interactions
print *, ""
print *, "Output files:"
print *, "  positions.dat - all particle positions"
print *, "  means.dat - mean square displacement"
print *, "  interactions.dat - detailed interaction log"
print *, "  traj_*.dat - individual trajectories"

! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,vhx,vhy,x0,y0,outside)

! Close files
close(11)
close(12)
close(13)
close(20)
close(21)
close(22)

end program main
