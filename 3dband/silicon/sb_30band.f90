subroutine band_30(h,kx,ky,kz,a)
implicit none 
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 

real(8)::El25,El2,E15,Eu1,El1,E12,Eu25,Eu2
real(8)::so25,so15
real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1
real(8),intent(in)::kx,ky,kz,h
real(8)::x,y,b,pi
complex(8)::zr,zi

El25=0.0d0
El2=4.185d0/13.6d0
E15=3.40d0/13.6d0
Eu1=7.07d0/13.6d0
El1=-12.92d0/13.6d0
E12=9.66d0/13.6d0
Eu25=12.78d0/13.6d0
Eu2=13.46d0/13.6d0
P=1.211d0
Q=1.044d0
R=0.574d0
P2=0.542d0
P1=0.296d0
Q1=0.742d0
R1=0.851d0
P3=1.236d0
T=1.097d0
T1=0.283d0

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
         a(i,i)=(El25+kx**2+ky**2+kz**2-so25)*zr
      else if (7<=i.and.i<=12) then
         a(i,i)=(E15+kx**2+ky**2+kz**2-so15)*zr
      else if (13<=i.and.i<=18) then
         a(i,i)=(Eu25+kx**2+ky**2+kz**2)*zr
      else if (19<=i.and.i<=20) then
         a(i,i)=(Eu1+kx**2+ky**2+kz**2)*zr
      else if (21<=i.and.i<=22) then
         a(i,i)=(El1+kx**2+ky**2+kz**2)*zr
      else if (23<=i.and.i<=24) then
         a(i,i)=(El2+kx**2+ky**2+kz**2)*zr
      else if (25<=i.and.i<=28) then
         a(i,i)=(E12+kx**2+ky**2+kz**2)*zr
      else        
         a(i,i)=(Eu2+kx**2+ky**2+kz**2)*zr
      end if
   end do
end do
!spin orbit matrix
   a(1,2)=so25*zi
   a(1,6)=so25*zr
   a(2,1)=-so25*zi
   a(2,6)=so25*zi
   a(3,4)=-so25*zr
   a(3,5)=-so25*zi
   a(4,5)=-so25*zi

   a(7,8)=so15*zi
   a(7,12)=so15*zr
   a(8,7)=-so15*zi
   a(8,12)=so15*zi
   a(9,10)=-so15*zr
   a(9,11)=-so15*zi
   a(10,11)=-so25*zi

   a(1,8)=Q*kz*zr
   a(1,9)=Q*ky*zr
   a(1,23)=P*kx*zr
   a(1,25)=-R*kx*zr
   a(1,26)=-sqrt(b)*R*kx*zr
   a(1,29)=P2*kx*zr
   
   a(2,7)=Q*kz*zr
   a(2,9)=Q*kx*zr
   a(2,23)=P*ky*zr
   a(2,25)=-R*ky*zr
   a(2,26)=sqrt(b)*R*ky*zr
   a(2,29)=P2*ky*zr

   a(3,7)=Q*ky*zr
   a(3,8)=Q*kx*zr
   a(3,23)=P*kz*zr
   a(3,25)=2.0d0*R*kz*zr
   a(3,29)=P2*kz*zr

   a(4,11)=Q*kz*zr
   a(4,12)=Q*ky*zr
   a(4,24)=P*kx*zr
   a(4,27)=-R*kx*zr
   a(4,28)=-sqrt(b)*R*kx*zr
   a(4,30)=P2*kx*zr

   a(5,10)=Q*kz*zr
   a(5,12)=Q*kx*zr
   a(5,24)=P*ky*zr
   a(5,27)=-R*ky*zr
   a(5,28)=sqrt(b)*R*ky*zr
   a(5,30)=P2*ky*zr

   a(6,10)=Q*ky*zr
   a(6,11)=Q*kx*zr
   a(6,24)=P*kz*zr
   a(6,27)=R*kz*2.0d0*zr
   a(6,30)=P2*kz*zr

   a(7,14)=-Q1*kz*zr
   a(7,15)=-Q1*ky*zr
   a(7,19)=T*kx*zr
   a(7,21)=T1*kx*zr

   a(8,13)=-Q1*kz*zr
   a(8,15)=-Q1*kx*zr
   a(8,19)=T*ky*zr
   a(8,21)=T1*ky*zr

   a(9,13)=-Q1*ky*zr
   a(9,14)=-Q1*kx*zr
   a(9,19)=T*kz*zr
   a(9,21)=T1*kz*zr

   a(10,17)=-Q1*kz*zr
   a(10,18)=-Q1*ky*zr
   a(10,20)=T*kx*zr
   a(10,22)=T1*kx*zr

   a(11,16)=-Q1*kz*zr
   a(11,18)=-Q1*kx*zr
   a(11,20)=T*ky*zr
   a(11,22)=T1*ky*zr

   a(12,16)=-Q1*ky*zr
   a(12,17)=-Q1*kx*zr
   a(12,20)=T*kz*zr
   a(12,22)=T1*kz*zr

   a(13,23)=-P1*kx*zr
   a(13,25)=-R1*kx*zr
   a(13,26)=-sqrt(b)*R1*kx*zr
   a(13,29)=P3*kx*zr
   
   a(14,23)=-P1*ky*zr
   a(14,25)=-R1*ky*zr
   a(14,26)=sqrt(b)*R1*ky*zr
   a(14,29)=P3*ky*zr

   a(15,23)=-P1*kz*zr
   a(15,25)=2.0d0*R1*kz*zr
   a(15,29)=P3*kz*zr
   
   a(16,24)=-P1*kx*zr
   a(16,27)=-R1*kx*zr
   a(16,28)=-sqrt(b)*R1*kx*zr
   a(16,30)=P3*kx*zr

   a(17,24)=-P1*ky*zr
   a(17,27)=-R1*ky*zr
   a(17,28)=sqrt(b)*R1*ky*zr
   a(17,30)=P3*ky*zr
   
   a(18,24)=-P1*kz*zr
   a(18,27)=R1*kz*2.0d0*zr
   a(18,30)=P3*kz*zr

   do i=1,n
      do j=1,n
         a(j,i)=conjg(a(i,j))
      end do
   end do

end subroutine band_30

