subroutine fullband(atomnum,laser_ev,nt,wc,cl25,c15,cu25,cu1,cl1,cl2,c12,cu2)
implicit none
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 
real(8)::kx,ky,kz,pi
real(8)::h,s,ca,cb,cc,cd
integer,intent(in)::atomnum,nt
real(8),intent(in)::laser_ev
real(8)::lattice,kmax,moji,laser
real(8),dimension(4*nt,n)::wc
real(8),dimension(4*nt,n)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2
real(8),dimension(20)::energypara,dipolepara

laser=laser_ev/13.6055d0

call paras(atomnum,lattice,energypara,dipolepara)
write(*,*) atomnum

pi=3.141592d0
lattice=5.43d0/0.529177d0  !格子定数/ボーア半径
kmax=2*pi/lattice

kx=0.5d0*kmax
ky=0.5d0*kmax
kz=0.5d0*kmax

do k=1,nt
   h=0.5d0*kmax*(k-1)/nt
   kx=0.5d0*kmax-h
   ky=0.5d0*kmax-h
   kz=0.5d0*kmax-h

call band_30(kx,ky,kz,energypara,dipolepara,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

   do i=1,n
      wc(k,i)=w(i)*13.6d0
   end do

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k,i)=cl25(k,i)+abs(a(j,i))**2
          write(90,*) j,abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k,i)=c15(k,i)+abs(a(j,i))**2
          write(90,*) j,abs(a(j,i))**2
         else if (13<=j.and.j<=18) then
           cu25(k,i)=cu25(k,i)+abs(a(j,i))**2
         else if (19<=j.and.j<=20) then
           cu1(k,i)=cu1(k,i)+abs(a(j,i))**2
         else if (21<=j.and.j<=22) then
           cl1(k,i)=cl1(k,i)+abs(a(j,i))**2
         else if (23<=j.and.j<=24) then
           cl2(k,i)=cl2(k,i)+abs(a(j,i))**2
         else if (25<=j.and.j<=28) then
           c12(k,i)=c12(k,i)+abs(a(j,i))**2
         else
           cu2(k,i)=cu2(k,i)+abs(a(j,i))**2
         end if
      end do
   end do

end do

kx=0     !初期化
ky=0
kz=0
a=0
do k=1,nt
   h=kmax*(k-1)/nt 
   kz=h
   call band_30(kx,ky,kz,energypara,dipolepara,a)
!call band_30(kx,ky,kz,a) 
call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) 
  do i=1,n
      wc(k+nt,i)=w(i)*13.6d0
  end do

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k+nt,i)=cl25(k+nt,i)+abs(a(j,i))**2
          write(90,*) j,abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k+nt,i)=c15(k+nt,i)+abs(a(j,i))**2
          write(90,*) j,abs(a(j,i))**2
         else if (13<=j.and.j<=18) then
           cu25(k+nt,i)=cu25(k+nt,i)+abs(a(j,i))**2
         else if (19<=j.and.j<=20) then
           cu1(k+nt,i)=cu1(k+nt,i)+abs(a(j,i))**2
         else if (21<=j.and.j<=22) then
           cl1(k+nt,i)=cl1(k+nt,i)+abs(a(j,i))**2
         else if (23<=j.and.j<=24) then
           cl2(k+nt,i)=cl2(k+nt,i)+abs(a(j,i))**2
         else if (25<=j.and.j<=28) then
           c12(k+nt,i)=c12(k+nt,i)+abs(a(j,i))**2
         else
           cu2(k+nt,i)=cu2(k+nt,i)+abs(a(j,i))**2
         end if
      end do
   end do
end do

end subroutine fullband
