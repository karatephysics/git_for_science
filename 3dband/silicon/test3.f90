subroutine ham(atomnum,nkx,nky,nkz,nt,energy,eigenv)
implicit none

character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n),work(lwork)
integer,intent(in)::atomnum,nkx,nky,nkz,nt
real(8)::kx,ky,kz,lattice,egap,pi
real(8),dimension(n)::para
real(8),dimension(n),intent(out)::energy
real(8),dimension(n,8),intent(out)::eigenv
real(8),dimension(20)::energypara,dipolepara

call paras(atomnum,lattice,energypara,dipolepara)

! aが固有関数なので注意
energy=0
pi=3.141592d0
eigenv=0

kx=(2*pi/lattice)*nkx/(nt-1)
ky=(2*pi/lattice)*nky/(nt-1)
kz=(2*pi/lattice)*nkz/(nt-1)

call band_30(kx,ky,kz,energypara,dipolepara,a)
call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
do i=1,n
    energy(i)=w(i)*13.6d0
end do
write(222,*) kx,ky,kz,energy(9)

do i=1,n
    do j=1,n
        if (1<=j.and.j<=6) then
        eigenv(i,1)=eigenv(i,1)+abs(a(j,i))**2
        else if (7<=j.and.j<=12) then
        eigenv(i,2)=eigenv(i,2)+abs(a(j,i))**2
        else if (13<=j.and.j<=18) then
        eigenv(i,3)=eigenv(i,3)+abs(a(j,i))**2
        else if (19<=j.and.j<=20) then
        eigenv(i,4)=eigenv(i,4)+abs(a(j,i))**2
        else if (21<=j.and.j<=22) then
        eigenv(i,5)=eigenv(i,5)+abs(a(j,i))**2
        else if (23<=j.and.j<=24) then
        eigenv(i,6)=eigenv(i,6)+abs(a(j,i))**2
        else if (25<=j.and.j<=28) then
        eigenv(i,7)=eigenv(i,7)+abs(a(j,i))**2
        else 
        eigenv(i,8)=eigenv(i,8)+abs(a(j,i))**2
        end if 
    end do
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
