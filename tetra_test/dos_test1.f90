subroutine parts_cal(e0,e1,e2,e3,VT,DOS,NOS)
implicit none
real(8),intent(in)::e0,e1,e2,e3,VT
integer::aa,i,j,j0,j1,j2,j3
real(8)::sort1,sort2
real(8)::e_ini,e_fin,nn,sa,VG
real(8)::e10,e20,e30,e21,e31,e32
real(8),dimension(4)::DOSpara
real(8),dimension(4)::stock
real(8),dimension(10001)::DOS,DOS_parts,NOS

!DOSpara(1)=51
DOSpara(1)=101
DOSpara(2)=10001
DOSpara(3)=0.7d0
DOSpara(4)=10.0d0
nn=DOSpara(2)

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

e_ini=DOSpara(3)
e_fin=DOSpara(4)
nn=DOSpara(2)
sa=(e_fin  -  e_ini)/(nn-1)

j0=nint((stock(1)-e_ini)/sa)
j1=nint((stock(2)-e_ini)/sa)
j2=nint((stock(3)-e_ini)/sa)
j3=nint((stock(4)-e_ini)/sa)

!write(23,*) j0,j1,j2,j3

VG=1/2.0d0/3.0d0*(DOSpara(1)-1)**3
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

!do i=1,10000
!   write(3,*) DOS(i)
!end do

end subroutine parts_cal
