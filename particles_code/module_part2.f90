module sectors
   use globals
   use Langevin
   implicit none
   integer :: M

contains
   !intakes coordinates and outputs index
   integer function coord_to_cell(px, py) result(idx)
      double precision, intent(in) :: px, py
      integer :: ix, iy

      M = floor(L / rc)

      ix = floor((px + L/2.0d0) / rc)
      iy = floor((py + L/2.0d0) / rc)

      !for particles on the boundary
      ix = min(max(ix, 0), M-1)
      iy = min(max(iy, 0), M-1)

      idx = iy * M + ix

   end function coord_to_cell

end module sectors
