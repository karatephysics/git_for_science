subroutine ham(a,b,c,energy)
implicit none

character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::aa(lda,n),work(lwork)

integer,intent(in)::a,b,c
real(8),dimension(30)::energy
real(8),dimension(2)::DOSpara
real(8)::kx,ky,kz,a0,egap,pi

DOSpara(1)=101

energy=0
pi=3.141592d0
a0=5.43d0/0.529d0

kx=(pi/a0)*a/(DOSpara(1)-1)
ky=(pi/a0)*b/(DOSpara(1)-1)
kz=(pi/a0)*c/(DOSpara(1)-1)

call band_30(kx,ky,kz,aa)
call zheev(jobz, uplo, n, aa, lda, w, work, lwork, rwork, info)
   do i=1,n
      energy(i)=w(i)*13.6d0
   end do

end subroutine ham

