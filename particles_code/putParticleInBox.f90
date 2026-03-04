   ! intakes coordinates and outputs index
   integer function indexfxn(px, py) result(idx)
     double precision, intent(in) :: px, py
     double precision :: boxSize,yCoord
      integer :: ix, iy

      ! the sector size      
      ! M = floor(L / rc)
      yCoord = - py
      ! coordinate + half length of boundary [-L/2,L/2] -> [0,L] 
      ! divided by effective particle size to give coordinate
      boxSize=L/real(b,8)
      ix = floor((px + L/2.0d0)/boxSize)
      iy = floor((yCoord + L/2.0d0)/boxSize)

      ! for particles on the boundary
      ix = min(max(ix, 0), b-1)
      iy = min(max(iy, 0), b-1)
      if (ix < 0 .or. ix >= b .or. iy < 0 .or. iy >= b) then
         print *, "ERROR: particle outside sector grid"
         print *, "px, py =", px, py
         stop
      end if

      ! call rey's function
      idx = TWOtoONE(ix,iy)

   end function indexfxn
