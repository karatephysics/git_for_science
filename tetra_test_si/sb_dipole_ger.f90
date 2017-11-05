subroutine dipole(nkx,nky,nkz,c)
implicit none
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k, l, m
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 
real(8)::kx,ky,kz
real(8)::h,s,ca,cb,cc,cd
integer,intent(in)::nkx,nky,nkz
!real(8),intent(in)::laser_ev
!integer,intent(in)::nt
real(8)::pi,a0,moji,laser,gap,energy
real(8),dimension(4)::DOSpara
real(8),dimension(10001)::c
real(8),dimension(n,n)::wc

real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1
real(8)::e_ini,e_fin,nn,sa
complex(8),dimension(100)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2
complex(8),dimension(10001)::e
complex(8)::zr,zi

P=1.345d0
Q=1.139d0
R=0.619d0
P2=0.429d0
P1=0.019d0
Q1=0.948d0
R1=1.076d0
P3=1.424d0
T=1.145d0
T1=0.281d0

zr=(1.0d0,0.0d0)
zi=(0.0d0,1.0d0)

DOSpara(1)=51!51
!DOSpara(1)=201
DOSpara(2)=10001
DOSpara(3)=0.8
DOSpara(4)=3
pi=3.141592d0
a0=5.6754d0/0.529177d0  !格子定数/ボーア半径

kx=(pi/a0)*nkx/(DOSpara(1)-1)
ky=(pi/a0)*nky/(DOSpara(1)-1)
kz=(pi/a0)*nkz/(DOSpara(1)-1)
w=0

call band_30(kx,ky,kz,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

cl25=0
c15=0
cu25=0
cu1=0
cl1=0
cl2=0
c12=0
cu2=0

e_ini=DOSpara(3)
e_fin=DOSpara(4)
nn=DOSpara(2)
sa=(e_fin-e_ini)/(nn-1)

  
do k=1,10000
energy=e_ini+sa*(k-1)
 do i=3,16
       do j=3,16
       gap=(w(i)-w(j))*13.6d0
       moji=gap-energy
       if (0<=moji.and.moji<=0.001) then
               do l=1,n
               do m=1,n
         if (1<=l.and.l<=6.and.23<=m.and.m<=24) then
         e(k)=e(k)+P*a(l,i)*a(m,j)
         e(k)=e(k)+P*a(l,j)*a(m,i)
         else if (13<=l.and.l<=18.and.23<=m.and.m<=24) then
         e(k)=e(k)+P1*a(l,i)*a(m,j)
         e(k)=e(k)+P1*a(l,j)*a(m,i)
         else if (1<=l.and.l<=6.and.29<=m.and.m<=30) then
         e(k)=e(k)+P2*a(l,i)*a(m,j)
         e(k)=e(k)+P2*a(l,j)*a(m,i)
         else if (13<=l.and.l<=18.and.29<=m.and.m<=30) then
         e(k)=e(k)+P3*a(l,i)*a(m,j)
         e(k)=e(k)+P3*a(l,j)*a(m,i)
         else if (1<=l.and.l<=6.and.7<=m.and.m<=12) then
         e(k)=e(k)+Q*a(l,i)*a(m,j)
         e(k)=e(k)+Q*a(l,j)*a(m,i)
         else if (13<=l.and.l<=18.and.7<=m.and.m<=12) then
         e(k)=e(k)+Q1*a(l,i)*a(m,j)
         e(k)=e(k)+Q1*a(l,j)*a(m,i)
         else if (1<=l.and.l<=6.and.25<=m.and.m<=28) then
         e(k)=e(k)+R*a(l,i)*a(m,j)
         e(k)=e(k)+R*a(l,j)*a(m,i)
         else if (13<=l.and.l<=18.and.25<=m.and.m<=28) then
         e(k)=e(k)+R1*a(l,i)*a(m,j)
         e(k)=e(k)+R1*a(l,j)*a(m,i)
         else if (21<=l.and.l<=22.and.7<=m.and.m<=12) then
         e(k)=e(k)+T*a(l,i)*a(m,j)
         e(k)=e(k)+T*a(l,j)*a(m,i)
         else if (19<=l.and.l<=20.and.7<=m.and.m<=12) then
         e(k)=e(k)+T1*a(l,i)*a(m,j)
         e(k)=e(k)+T1*a(l,j)*a(m,i)
         end if
       end do
       end do
       end if
       end do
  end do
end do

do k=1,10000
        c(k)=abs(e(k))**2
end do

end subroutine dipole 
