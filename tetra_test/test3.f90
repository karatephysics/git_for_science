subroutine parts_cal(para,e0,e1,e2,e3,VT,DOS,NOS)
implicit none
real(8),intent(in)::e0,e1,e2,e3,VT
real(8),dimension(4),intent(in)::para
integer::aa,i,j,j0,j1,j2,j3
real(8)::sort1,sort2
real(8)::e_ini,e_fin,nn,sa,VG
real(8)::e10,e20,e30,e21,e31,e32
real(8),dimension(4)::stock
real(8),dimension(2001)::DOS,DOS_parts,NOS

stock(1)=e0
stock(2)=e1
stock(3)=e2
stock(4)=e3

do i=1,4
   do j=i+1,4
      if (stock(i)>stock(j)) then
         sort1=stock(i)
         sort2=stock(j)
         stock(i)=sort2
         stock(j)=sort1
      end if
   end do
end do

!do i=1,4
!   write(14,*) i,stock(i)
!end do

nn=para(2)
e_ini=para(3)
e_fin=para(4)
sa=(e_fin  -  e_ini)/(nn-1)

j0=nint((stock(1)-e_ini)/sa)
j1=nint((stock(2)-e_ini)/sa)
j2=nint((stock(3)-e_ini)/sa)
j3=nint((stock(4)-e_ini)/sa)

!write(23,*) j0,j1,j2,j3

VG=1/2.0d0/3.0d0*(para(1)-1)**3
e10=stock(2)-stock(1)
e20=stock(3)-stock(1)
e30=stock(4)-stock(1)
e21=stock(3)-stock(2)
e31=stock(4)-stock(2)
e32=stock(4)-stock(3)
!e10=e1-e0
!e20=e2-e0
!e30=e3-e0
!e21=e2-e1
!e31=e3-e1
!e32=e3-e2
aa=0

!write(24,'(6f9.5)') e10,e20,e30,e21,e31,e32

if (stock(1)==0.and.stock(2)==0.and.stock(3)==0.and.stock(4)==0) then
else
   do aa=j0+1,j1
      if(aa<0.or.aa>(nn-1)) then
      else
        DOS(aa)=DOS(aa)+VT/VG*3*(aa*sa+e_ini-stock(1))**2/(e10*e20*e30)
         NOS(aa)=NOS(aa)+VT/VG*(aa*sa+e_ini  -  stock(1))**3/(e10*e20*e30)
!         DOS_parts(aa)=DOS_parts(aa)+VT/VG*3*(aa*sa+e_ini  -  e0)**2/(e10*e20*e30)
!      write(25,*) aa,DOS(aa)
!     write(26,*) VT,VG,sa
!      write(26,*) e_ini,e0,e10,e20,e30

      end if
   end do

   do aa=j1+1,j2
      if (aa<0.or.aa>(nn-1)) then
      else
         DOS(aa)=DOS(aa)+VT/VG/(e20*e30)*(3*e10 +6*(aa*sa+e_ini-stock(2))&
         -3*(e20+e31)/(e21*e31)*(aa*sa+e_ini-stock(2))**2 )
         NOS(aa)=NOS(aa)+VT/VG*1/(e20*e30)*(e10**2+3*e10*(aa*sa+e_ini-stock(3))&
         +3*(aa*sa+e_ini-stock(3))**2-(e20 + e31)/(e21*e31)*(aa*sa+e_ini-stock(3))**3 )
!         DOS_parts(aa)=DOS_part(aa)+VT/VG/(e20*e30)*(3*e10 +6*(aa*sa+e_ini-e1)&
!         -3*(e20+e31)/(e21*e31)*(aa*sa+e_ini-e1)**2 )
      end if
   end do

   do aa=j2+1,j3
      if(aa<0.or.aa>(nn-1)) then
      else
         DOS(aa)=DOS(aa)+VT/VG*( 3*(stock(4)  -  (aa*sa+e_ini))**2/(e30*e31*e32) )
         NOS(aa)=NOS(aa)+VT*VG*(1  -(stock(4) -  (aa*sa+e_ini))**3/(e30*e31*e32)  )
!         DOS_parts(aa)=DOS_parts(aa)+VT/VG*( 3*(e3  -  (aa*sa+e_ini))**2/(e30*e31*e32) )
      end if
   end do
end if

   write(50,'(2f15.10)') DOS(1050),NOS(1050)
   write(51,'(2f15.10)') DOS(1051),NOS(1051)
   write(52,'(2f15.10)') DOS(1052),NOS(1052)
   write(53,'(2f15.10)') DOS(1053),NOS(1053)
   write(54,'(2f15.10)') DOS(1054),NOS(1054)
   write(55,'(2f15.10)') DOS(1055),NOS(1055)
   write(56,'(2f15.10)') DOS(1056),NOS(1056)

   write(57,'(2f15.10)') DOS(257),NOS(257)
   write(58,'(2f15.10)') DOS(258),NOS(258)
   write(59,'(2f15.10)') DOS(259),NOS(259)
   write(60,'(2f15.10)') DOS(260),NOS(260)
   write(61,'(2f15.10)') DOS(261),NOS(261)
   write(62,'(2f15.10)') DOS(262),NOS(262)
   write(63,'(2f15.10)') DOS(263),NOS(263)
   write(64,'(2f15.10)') DOS(264),NOS(264)

end subroutine parts_cal
