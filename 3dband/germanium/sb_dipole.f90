subroutine ham(a,b,c,nt,energy,eigenv)
implicit none

character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::aa(lda,n),work(lwork)
integer,intent(in)::a,b,c,nt
real(8)::kx,ky,kz,a0,egap,pi
real(8),dimension(n)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2
real(8),dimension(n)::para
real(8),dimension(n),intent(out)::energy
real(8),dimension(n,8),intent(out)::eigenv

! aaが固有関数なので注意
energy=0
pi=3.141592d0
a0=5.6754d0/0.529d0
cl25=0
cu25=0
c15=0
cu1=0
cl1=0
cl2=0
c12=0
cu2=0 

kx=(pi/a0)*a/(nt-1)
ky=(pi/a0)*b/(nt-1)
kz=(pi/a0)*c/(nt-1)

call band_30(kx,ky,kz,aa)
call zheev(jobz, uplo, n, aa, lda, w, work, lwork, rwork, info)
do i=1,n
    energy(i)=w(i)*13.6d0
end do

do i=1,n
    do j=1,n
        if (1<=j.and.j<=6) then
        cl25(i)=cl25(i)+abs(aa(j,i))**2
        else if (7<=j.and.j<=12) then
        c15(i)=c15(i)+abs(aa(j,i))**2
        else if (13<=j.and.j<=18) then
        cu25(i)=cu25(i)+abs(aa(j,i))**2
        else if (19<=j.and.j<=20) then
        cu1(i)=cu1(i)+abs(aa(j,i))**2
        else if (21<=j.and.j<=22) then
        cl1(i)=cl1(i)+abs(aa(j,i))**2
        else if (23<=j.and.j<=24) then
        cl2(i)=cl2(i)+abs(aa(j,i))**2
        else if (25<=j.and.j<=28) then
        c12(i)=c12(i)+abs(aa(j,i))**2
        else 
        cu2(i)=cu2(i)+abs(aa(j,i))**2
        end if 
    end do
end do

do i=1,n
    eigenv(i,1)=cl25(i)
    eigenv(i,2)=c15(i)
    eigenv(i,3)=cu25(i)
    eigenv(i,4)=cu1(i)
    eigenv(i,5)=cl1(i)
    eigenv(i,6)=cl2(i)
    eigenv(i,7)=c12(i)
    eigenv(i,8)=cu2(i)
end do

!規格化条件のチェック
do i=1,n
    para=0
    do j=1,8
        para(i)=para(i)+eigenv(i,j)
    end do
    if (nint(para(i))==1) then
        write(202,*) para(i)
        else
        write(*,*) para(i)
    end if 
end do

end subroutine ham
