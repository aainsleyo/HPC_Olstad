module sectors
   use globals
   use Langevin
   use griding ! Rey's TWOtoONE
   implicit none

contains
   ! intakes coordinates and outputs index
   integer function indexfxn(px, py) result(idx)
      double precision, intent(in) :: px, py
      integer :: ix, iy

      ! the sector size      
      ! M = floor(L / rc)

      ! coordinate + half length of boundary [-L/2,L/2] -> [0,L] 
      ! divided by effective particle size to give coordinate
      ix = floor((px + L/2.0d0) / rc)
      iy = floor((py + L/2.0d0) / rc)

      ! for particles on the boundary
      !ix = min(max(ix, 0), M-1)
      !iy = min(max(iy, 0), M-1)
      if (ix < 0 .or. ix >= M .or. iy < 0 .or. iy >= M) then
         print *, "ERROR: particle outside sector grid"
         print *, "px, py =", px, py
         stop
      end if

      ! call rey's function
      idx = TWOtoONE(ix,iy,M)

   end function indexfxn

end module sectors
