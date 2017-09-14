subroutine band_30(kx,ky,kz,energypara,dipolepara,a)
implicit none 
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 

real(8),intent(in)::kx,ky,kz
real(8),dimension(20),intent(in)::energypara
real(8),dimension(20),intent(in)::dipolepara
real(8)::b,pi
complex(8)::zr,zi

b=3.0d0
pi=3.141592d0 
zr=(1.0d0,0.0d0) 
zi=(0.0d0,1.0d0)

do i=1,n
   do j=1,n
      a(i,j)=0.0d0
   end do
end do

do i=1,n
   do j=1,n
      if (1<=i.and.i<=6) then
         a(i,i)=(energypara(2)+kx**2+ky**2+kz**2-energypara(11))*zr
      else if (7<=i.and.i<=12) then
         a(i,i)=(energypara(4)+kx**2+ky**2+kz**2-energypara(12))*zr
      else if (13<=i.and.i<=18) then
         a(i,i)=(energypara(7)+kx**2+ky**2+kz**2)*zr
      else if (19<=i.and.i<=20) then
         a(i,i)=(energypara(5)+kx**2+ky**2+kz**2)*zr
      else if (21<=i.and.i<=22) then
         a(i,i)=(energypara(1)+kx**2+ky**2+kz**2)*zr
      else if (23<=i.and.i<=24) then
         a(i,i)=(energypara(3)+kx**2+ky**2+kz**2)*zr
      else if (25<=i.and.i<=28) then
         a(i,i)=(energypara(6)+kx**2+ky**2+kz**2)*zr
      else        
         a(i,i)=(energypara(8)+kx**2+ky**2+kz**2)*zr
      end if
   end do
end do
!spin orbit matrix
   a(1,2)=energypara(11)*zi
   a(1,6)=energypara(11)*zr
   a(2,1)=-energypara(11)*zi
   a(2,6)=energypara(11)*zi
   a(3,4)=-energypara(11)*zr
   a(3,5)=-energypara(11)*zi
   a(4,5)=-energypara(11)*zi

   a(7,8)=energypara(12)*zi
   a(7,12)=energypara(12)*zr
   a(8,7)=-energypara(12)*zi
   a(8,12)=energypara(12)*zi
   a(9,10)=-energypara(12)*zr
   a(9,11)=-energypara(12)*zi
   a(10,11)=-energypara(11)*zi

   a(1,8)=dipolepara(5)*kz*zr
   a(1,9)=dipolepara(5)*ky*zr
   a(1,23)=dipolepara(1)*kx*zr
   a(1,25)=-dipolepara(7)*kx*zr
   a(1,26)=-sqrt(b)*dipolepara(7)*kx*zr
   a(1,29)=dipolepara(3)*kx*zr
   
   a(2,7)=dipolepara(5)*kz*zr
   a(2,9)=dipolepara(5)*kx*zr
   a(2,23)=dipolepara(1)*ky*zr
   a(2,25)=-dipolepara(7)*ky*zr
   a(2,26)=sqrt(b)*dipolepara(7)*ky*zr
   a(2,29)=dipolepara(3)*ky*zr

   a(3,7)=dipolepara(5)*ky*zr
   a(3,8)=dipolepara(5)*kx*zr
   a(3,23)=dipolepara(1)*kz*zr
   a(3,25)=2.0d0*dipolepara(7)*kz*zr
   a(3,29)=dipolepara(3)*kz*zr

   a(4,11)=dipolepara(5)*kz*zr
   a(4,12)=dipolepara(5)*ky*zr
   a(4,24)=dipolepara(1)*kx*zr
   a(4,27)=-dipolepara(7)*kx*zr
   a(4,28)=-sqrt(b)*dipolepara(7)*kx*zr
   a(4,30)=dipolepara(3)*kx*zr

   a(5,10)=dipolepara(5)*kz*zr
   a(5,12)=dipolepara(5)*kx*zr
   a(5,24)=dipolepara(1)*ky*zr
   a(5,27)=-dipolepara(7)*ky*zr
   a(5,28)=sqrt(b)*dipolepara(7)*ky*zr
   a(5,30)=dipolepara(3)*ky*zr

   a(6,10)=dipolepara(5)*ky*zr
   a(6,11)=dipolepara(5)*kx*zr
   a(6,24)=dipolepara(1)*kz*zr
   a(6,27)=dipolepara(7)*kz*2.0d0*zr
   a(6,30)=dipolepara(3)*kz*zr

   a(7,14)=-dipolepara(6)*kz*zr
   a(7,15)=-dipolepara(6)*ky*zr
   a(7,19)=dipolepara(9)*kx*zr
   a(7,21)=dipolepara(10)*kx*zr

   a(8,13)=-dipolepara(6)*kz*zr
   a(8,15)=-dipolepara(6)*kx*zr
   a(8,19)=dipolepara(9)*ky*zr
   a(8,21)=dipolepara(10)*ky*zr

   a(9,13)=-dipolepara(6)*ky*zr
   a(9,14)=-dipolepara(6)*kx*zr
   a(9,19)=dipolepara(9)*kz*zr
   a(9,21)=dipolepara(10)*kz*zr

   a(10,17)=-dipolepara(6)*kz*zr
   a(10,18)=-dipolepara(6)*ky*zr
   a(10,20)=dipolepara(9)*kx*zr
   a(10,22)=dipolepara(10)*kx*zr

   a(11,16)=-dipolepara(6)*kz*zr
   a(11,18)=-dipolepara(6)*kx*zr
   a(11,20)=dipolepara(9)*ky*zr
   a(11,22)=dipolepara(10)*ky*zr

   a(12,16)=-dipolepara(6)*ky*zr
   a(12,17)=-dipolepara(6)*kx*zr
   a(12,20)=dipolepara(9)*kz*zr
   a(12,22)=dipolepara(10)*kz*zr

   a(13,23)=-dipolepara(2)*kx*zr
   a(13,25)=-dipolepara(8)*kx*zr
   a(13,26)=-sqrt(b)*dipolepara(8)*kx*zr
   a(13,29)=dipolepara(4)*kx*zr
   
   a(14,23)=-dipolepara(2)*ky*zr
   a(14,25)=-dipolepara(8)*ky*zr
   a(14,26)=sqrt(b)*dipolepara(8)*ky*zr
   a(14,29)=dipolepara(4)*ky*zr

   a(15,23)=-dipolepara(2)*kz*zr
   a(15,25)=2.0d0*dipolepara(8)*kz*zr
   a(15,29)=dipolepara(4)*kz*zr
   
   a(16,24)=-dipolepara(2)*kx*zr
   a(16,27)=-dipolepara(8)*kx*zr
   a(16,28)=-sqrt(b)*dipolepara(8)*kx*zr
   a(16,30)=dipolepara(4)*kx*zr

   a(17,24)=-dipolepara(2)*ky*zr
   a(17,27)=-dipolepara(8)*ky*zr
   a(17,28)=sqrt(b)*dipolepara(8)*ky*zr
   a(17,30)=dipolepara(4)*ky*zr
   
   a(18,24)=-dipolepara(2)*kz*zr
   a(18,27)=dipolepara(8)*kz*2.0d0*zr
   a(18,30)=dipolepara(4)*kz*zr

   do i=1,n
      do j=1,n
         a(j,i)=conjg(a(i,j))
      end do
   end do

end subroutine band_30

