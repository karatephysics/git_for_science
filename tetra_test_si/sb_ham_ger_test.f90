subroutine ham(nkx,nky,nkz,wc)
implicit none
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 
real(8)::kx,ky,kz
real(8)::h,s,ca,cb,cc,cd
real(8),intent(in)::nkx,nky,nkz
!real(8),intent(in)::laser_ev
!integer,intent(in)::nt
real(8)::pi,a0,moji,laser
real(8),dimension(4)::DOSpara
real(8),dimension(100)::wc

DOSpara(1)=51
!DOSpara(1)=101
pi=3.141592d0
a0=5.6754d0/0.529177d0  !格子定数/ボーア半径

kx=(pi/a0)*nkx/(DOSpara(1)-1)
ky=(pi/a0)*nky/(DOSpara(1)-1)
kz=(pi/a0)*nkz/(DOSpara(1)-1)
w=0
wc=0

call band_30(kx,ky,kz,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

k=1
   do i=9,20
        do j=5,8
        wc(k)=w(i)-w(j)
        k=k+1
        end do
   end do
cl25=0
c15=0
cu25=0
cu1=0
cl1=0
cl2=0
c12=0
   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(i)=cl25(i)+abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(i)=c15(i)+abs(a(j,i))**2
         else if (13<=j.and.j<=18) then
           cu25(i)=cu25(i)+abs(a(j,i))**2
         else if (19<=j.and.j<=20) then
           cu1(i)=cu1(i)+abs(a(j,i))**2
         else if (21<=j.and.j<=22) then
           cl1(i)=cl1(i)+abs(a(j,i))**2
         else if (23<=j.and.j<=24) then
           cl2(i)=cl2(i)+abs(a(j,i))**2
         else if (25<=j.and.j<=28) then
           c12(i)=c12(i)+abs(a(j,i))**2
        else
           cu2(i)=cu2(i)+abs(a(j,i))**2
         end if
      end do
   end do

   do i=1,n
       ! write(212,"(2f11.2)") real(cl25(i)),real(cu25(i)
)
        write(212,"(8f11.4)") real(cl25(i)),real(cu25(i))
,real(c15(i)),real(cu1(i)),real(cl1(i))&
                ,real(cl2(i)),real(c12(i)),real(cu2(i))

!   write(*,*) k
end subroutine ham
