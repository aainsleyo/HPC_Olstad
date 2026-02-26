module globals
implicit none
double precision :: L = 10.0d0
double precision :: rc = 2.0d0
integer :: M
contains
subroutine update_M()
    M = floor(L/rc)
end subroutine update_M
end module globals

module griding
implicit none
contains
function TWOtoONE(x,y,box_size) result(id)
    integer,intent(in) :: x,y,box_size
    integer :: id
    id = y*box_size + x
end function TWOtoONE
end module griding

module sectors
use globals
use griding
implicit none
contains

integer function indexfxn(px,py) result(idx)
    double precision,intent(in) :: px,py
    integer :: ix,iy

    ix = floor((px + L/2d0)/rc)
    iy = floor((py + L/2d0)/rc)

    if (ix < 0 .or. ix >= M .or. iy < 0 .or. iy >= M) then
        print *, "ERROR: Out of bounds"
        stop
    end if

    idx = TWOtoONE(ix,iy,M)
end function indexfxn

end module sectors

program test_particle
use globals
use sectors
implicit none

double precision :: px, py
integer :: idx

call update_M()

! TEST PARTICLE
px = 1.3d0
py = -2.7d0

idx = indexfxn(px,py)

print *, "Testing particle:"
print *, "px =", px
print *, "py =", py
print *, "Sector index =", idx

end program test_particle
