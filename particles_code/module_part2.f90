module sectors
   use globals
   use Langevin
   implicit none
   integer :: M

contains
   ! intakes coordinates and outputs index
   integer function coord_to_cell(px, py) result(idx)
      double precision, intent(in) :: px, py
      integer :: ix, iy

      ! the sector size      
      M = floor(L / rc)

      ! coordinate + half length of box divided by effective particle size
      ix = floor((px + L/2.0d0) / rc)
      iy = floor((py + L/2.0d0) / rc)

      ! for particles on the boundary
      ix = min(max(ix, 0), M-1)
      iy = min(max(iy, 0), M-1)

      ! standard 2D-to-1D mapping for typewriter order (row major order)
      idx = iy * M + ix

   end function coord_to_cell

end module sectors
