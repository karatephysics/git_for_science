program test

implicit none

integer::i,j
real::b1,b2
real,dimension(4)::a,b

a=0
b=0

a(1)=1
a(2)=3
a(3)=0
a(4)=5

do i=1,4
   do j=i+1,4
      if (a(i)>a(j)) then
         b1=a(i)
         b2=a(j)
         a(i)=b2
         a(j)=b1
      end if
   end do
end do

do i=1,4
   write(*,*) a(i)
end do

end program test
