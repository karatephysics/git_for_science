program test
implicit none
real(8)::e0,e1,e2,e3,VT
integer::aa,i,j,j0,j1,j2,j3
real(8)::sort1,sort2
real(8)::e_ini,e_fin,nn,sa,VG
real(8)::e10,e20,e30,e21,e31,e32
real(8),dimension(4)::stock

stock(1)=1!e0
stock(2)=3!e1
stock(3)=2!e2
stock(4)=5!e3

do i=1,4
   do j=i+1,4
      if (stock(i)>stock(j)) then
         sort1=stock(i)
         sort2=stock(j)
         stock(i)=sort1
         stock(j)=sort2
      end if
   end do
end do

do i=1,4
   write(*,*) i,stock(i)
end do

end program test
