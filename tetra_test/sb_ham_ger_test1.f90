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

!DOSpara(1)=51
DOSpara(1)=101
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
   do i=9,30,2
        do j=1,8,2
        wc(k)=w(i)-w(j)
        k=k+1
        end do
   end do

!   write(*,*) k
end subroutine ham
