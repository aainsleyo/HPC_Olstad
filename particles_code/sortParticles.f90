subroutine makeNR(nr,in)
  integer :: particle,sector,nr(0:b*b),in(n)
  nr = 0
  do particle=1,n
     sector=indexfxn(x(particle),y(particle))
     in(particle)=sector
     nr(sector)=nr(sector)+1
  end do
  return
end subroutine makeNR

subroutine order(x,y,vx,vy,x0,y0,lim)
! Order the list of positions by sector and find starting and final index for each sector
! In: x and y coordinates and velocities. Out: ordered lists x, y, vx and vy and array lim with one row 
! for each sector, first column is start index, second is end index so that particles lim(k,1) through lim(k,2) reside in sector k.
integer :: lim(0:b*b,2),in(n),nr(0:b*b),ct(0:b*b),k
double precision :: x(n),y(n),d1(n),d2(n),d3(n),d4(n),d5(n),d6(n),vx(n),vy(n),x0(n),y0(n)

call makeNR(nr,in)

! Set loop limits based on the number of particles in each sector
lim(0,1)=1
lim(0,2)=nr(0)

do k=1,b*b
   lim(k,1)=lim(k-1,2)+1
   lim(k,2)=lim(k-1,2)+nr(k)
end do

! Re-order particle list
d1=x
d2=y
d3=vx
d4=vy
d5=x0
d6=y0
ct=0
do k=1,n
   x(lim(in(k),1)+ct(in(k)))=d1(k)
   y(lim(in(k),1)+ct(in(k)))=d2(k)
   vx(lim(in(k),1)+ct(in(k)))=d3(k)
   vy(lim(in(k),1)+ct(in(k)))=d4(k)
   x0(lim(in(k),1)+ct(in(k)))=d5(k)
   y0(lim(in(k),1)+ct(in(k)))=d6(k)
   ct(in(k))=ct(in(k))+1
end do

end subroutine order
