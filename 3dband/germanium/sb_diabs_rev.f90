subroutine fullband(laser_ev,nt,wc,cl25,c15,cu25,cu1,cl1,cl2,c12,cu2)
implicit none
character, parameter :: jobz="v", uplo="l" 
integer, parameter ::n=30, lda=n, lwork=2*n-1 
integer ::info, i, j, k 
double precision, dimension(3*n-2) :: rwork 
double precision, dimension(n) ::w 
complex(kind(0d0)) ::a(lda,n), work(lwork) 
real(8)::kx,ky,kz,pi
real(8)::h,s,ca,cb,cc,cd
real(8),intent(in)::laser_ev
integer,intent(in)::nt
real(8)::ap,age,kmax,moji,laser
real(8),dimension(4*nt,n,n)::gap
real(8),dimension(4*nt,n)::wc
real(8),dimension(4*nt,n)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2

pi=3.141592d0
age=5.6754d0/0.529177d0  !格子定数/ボーア半径
ap=age
kmax=2*pi/ap
laser=laser_ev/13.6055d0
gap=0
kx=0.5d0*kmax
ky=0.5d0*kmax
kz=0.5d0*kmax

do k=1,nt
   h=0.5d0*kmax*(k-1)/nt
   kx=0.5d0*kmax-h
   ky=0.5d0*kmax-h
   kz=0.5d0*kmax-h
call band_30(kx,ky,kz,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

   do i=1,n
      wc(k,i)=w(i)*13.6d0
   end do

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k,i)=cl25(k,i)+abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k,i)=c15(k,i)+abs(a(j,i))**2
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
         else if (29<=j.and.j<=30) then
           cu2(k,i)=cu2(k,i)+abs(a(j,i))**2
         else
             write(*,*)"none"
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
call band_30(kx,ky,kz,a) 
call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) 
  do i=1,n
      wc(k+nt,i)=w(i)*13.6d0
  end do

if (211<=k) then
    wc(k+nt,9)=w(11)*13.6d0
    wc(k+nt,10)=w(12)*13.6d0
    wc(k+nt,11)=w(9)*13.6d0
    wc(k+nt,12)=w(10)*13.6d0

    ca=cl25(k+nt,11)
    cb=cl25(k+nt,12)
    cc=cl25(k+nt,9)
    cd=cl25(k+nt,10)
    cl25(k+nt,9)=ca
    cl25(k+nt,10)=cb
    cl25(k+nt,11)=cc
    cl25(k+nt,12)=cd

    ca=c15(k+nt,11)
    cb=c15(k+nt,12)
    cc=c15(k+nt,9)
    cd=c15(k+nt,10)
    c15(k+nt,9)=ca
    c15(k+nt,10)=cb
    c15(k+nt,11)=cc
    c15(k+nt,12)=cd

    ca=cu25(k+nt,11)
    cb=cu25(k+nt,12)
    cc=cu25(k+nt,9)
    cd=cu25(k+nt,10)
    cu25(k+nt,9)=ca
    cu25(k+nt,10)=cb
    cu25(k+nt,11)=cc
    cu25(k+nt,12)=cd

    ca=cu1(k+nt,11)
    cb=cu1(k+nt,12)
    cc=cu1(k+nt,9)
    cd=cu1(k+nt,10)
    cu1(k+nt,9)=ca
    cu1(k+nt,10)=cb
    cu1(k+nt,11)=cc
    cu1(k+nt,12)=cd

    ca=cl1(k+nt,11)
    cb=cl1(k+nt,12)
    cc=cl1(k+nt,9)
    cd=cl1(k+nt,10)
    cl1(k+nt,9)=ca
    cl1(k+nt,10)=cb
    cl1(k+nt,11)=cc
    cl1(k+nt,12)=cd

    ca=cl2(k+nt,11)
    cb=cl2(k+nt,12)
    cc=cl2(k+nt,9)
    cd=cl2(k+nt,10)
    cl2(k+nt,9)=ca
    cl2(k+nt,10)=cb
    cl2(k+nt,11)=cc
    cl2(k+nt,12)=cd

    ca=c12(k+nt,11)
    cb=c12(k+nt,12)
    cc=c12(k+nt,9)
    cd=c12(k+nt,10)
    c12(k+nt,9)=ca
    c12(k+nt,10)=cb
    c12(k+nt,11)=cc
    c12(k+nt,12)=cd

    ca=cu2(k+nt,11)
    cb=cu2(k+nt,12)
    cc=cu2(k+nt,9)
    cd=cu2(k+nt,10)
    cu2(k+nt,9)=ca
    cu2(k+nt,10)=cb
    cu2(k+nt,11)=cc
    cu2(k+nt,12)=cd

end if

!if (13<=k.or.485<=k) then
if (k<=485) then
    wc(k+nt,19)=w(21)*13.6d0
    wc(k+nt,20)=w(22)*13.6d0
    wc(k+nt,21)=w(19)*13.6d0
    wc(k+nt,22)=w(20)*13.6d0

    ca=cl25(k+nt,21)
    cb=cl25(k+nt,22)
    cc=cl25(k+nt,19)
    cd=cl25(k+nt,20)
    cl25(k+nt,21)=ca
    cl25(k+nt,22)=cb
    cl25(k+nt,19)=cc
    cl25(k+nt,20)=cd

    ca=c15(k+nt,21)
    cb=c15(k+nt,22)
    cc=c15(k+nt,19)
    cd=c15(k+nt,20)
    c15(k+nt,21)=ca
    c15(k+nt,22)=cb
    c15(k+nt,19)=cc
    c15(k+nt,20)=cd

    ca=cu25(k+nt,21)
    cb=cu25(k+nt,22)
    cc=cu25(k+nt,19)
    cd=cu25(k+nt,20)
    cu25(k+nt,21)=ca
    cu25(k+nt,22)=cb
    cu25(k+nt,19)=cc
    cu25(k+nt,20)=cd

    ca=cu1(k+nt,21)
    cb=cu1(k+nt,22)
    cc=cu1(k+nt,19)
    cd=cu1(k+nt,20)
    cu1(k+nt,21)=ca
    cu1(k+nt,22)=cb
    cu1(k+nt,19)=cc
    cu1(k+nt,20)=cd

    ca=cl1(k+nt,21)
    cb=cl1(k+nt,22)
    cc=cl1(k+nt,19)
    cd=cl1(k+nt,20)
    cl1(k+nt,21)=ca
    cl1(k+nt,22)=cb
    cl1(k+nt,19)=cc
    cl1(k+nt,20)=cd

    ca=cl2(k+nt,21)
    cb=cl2(k+nt,22)
    cc=cl2(k+nt,19)
    cd=cl2(k+nt,20)
    cl2(k+nt,21)=ca
    cl2(k+nt,22)=cb
    cl2(k+nt,19)=cc
    cl2(k+nt,20)=cd

    ca=c12(k+nt,21)
    cb=c12(k+nt,22)
    cc=c12(k+nt,19)
    cd=c12(k+nt,20)
    c12(k+nt,21)=ca
    c12(k+nt,22)=cb
    c12(k+nt,19)=cc
    c12(k+nt,20)=cd

    ca=cu2(k+nt,21)
    cb=cu2(k+nt,22)
    cc=cu2(k+nt,19)
    cd=cu2(k+nt,20)
    cu2(k+nt,21)=ca
    cu2(k+nt,22)=cb
    cu2(k+nt,19)=cc
    cu2(k+nt,20)=cd

end if

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k+nt,i)=cl25(k+nt,i)+abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k+nt,i)=c15(k+nt,i)+abs(a(j,i))**2
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

a=0
kx=kmax
ky=0.0d0
kz=0.0d0

do k=1,nt
   h=0.75*kmax*(k-1)/nt
   kx=kmax-h/3
   ky=h
call band_30(kx,ky,kz,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

  do i=1,n
     wc(k+nt*2,i)=w(i)*13.6d0
  end do

if (600<=k) then
    wc(k+nt*2,9)=w(11)*13.6d0
    wc(k+nt*2,10)=w(12)*13.6d0
    wc(k+nt*2,11)=w(9)*13.6d0
    wc(k+nt*2,12)=w(10)*13.6d0

    ca=cl25(k+nt*2,11)
    cb=cl25(k+nt*2,12)
    cc=cl25(k+nt*2,9)
    cd=cl25(k+nt*2,10)
    cl25(k+nt*2,9)=ca
    cl25(k+nt*2,10)=cb
    cl25(k+nt*2,11)=cc
    cl25(k+nt*2,12)=cd

    ca=c15(k+nt*2,11)
    cb=c15(k+nt*2,12)
    cc=c15(k+nt*2,9)
    cd=c15(k+nt*2,10)
    c15(k+nt*2,9)=ca
    c15(k+nt*2,10)=cb
    c15(k+nt*2,11)=cc
    c15(k+nt*2,12)=cd

    ca=cu25(k+nt*2,11)
    cb=cu25(k+nt*2,12)
    cc=cu25(k+nt*2,9)
    cd=cu25(k+nt*2,10)
    cu25(k+nt*2,9)=ca
    cu25(k+nt*2,10)=cb
    cu25(k+nt*2,11)=cc
    cu25(k+nt*2,12)=cd

    ca=cu1(k+nt*2,11)
    cb=cu1(k+nt*2,12)
    cc=cu1(k+nt*2,9)
    cd=cu1(k+nt*2,10)
    cu1(k+nt*2,9)=ca
    cu1(k+nt*2,10)=cb
    cu1(k+nt*2,11)=cc
    cu1(k+nt*2,12)=cd

    ca=cl1(k+nt*2,11)
    cb=cl1(k+nt*2,12)
    cc=cl1(k+nt*2,9)
    cd=cl1(k+nt*2,10)
    cl1(k+nt*2,9)=ca
    cl1(k+nt*2,10)=cb
    cl1(k+nt*2,11)=cc
    cl1(k+nt*2,12)=cd

    ca=cl2(k+nt*2,11)
    cb=cl2(k+nt*2,12)
    cc=cl2(k+nt*2,9)
    cd=cl2(k+nt*2,10)
    cl2(k+nt*2,9)=ca
    cl2(k+nt*2,10)=cb
    cl2(k+nt*2,11)=cc
    cl2(k+nt*2,12)=cd

    ca=c12(k+nt*2,11)
    cb=c12(k+nt*2,12)
    cc=c12(k+nt*2,9)
    cd=c12(k+nt*2,10)
    c12(k+nt*2,9)=ca
    c12(k+nt*2,10)=cb
    c12(k+nt*2,11)=cc
    c12(k+nt*2,12)=cd

    ca=cu2(k+nt*2,11)
    cb=cu2(k+nt*2,12)
    cc=cu2(k+nt*2,9)
    cd=cu2(k+nt*2,10)
    cu2(k+nt*2,9)=ca
    cu2(k+nt*2,10)=cb
    cu2(k+nt*2,11)=cc
    cu2(k+nt*2,12)=cd
end if

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k+nt*2,i)=cl25(k+nt*2,i)+abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k+nt*2,i)=c15(k+nt*2,i)+abs(a(j,i))**2
         else if (13<=j.and.j<=18) then
           cu25(k+nt*2,i)=cu25(k+nt*2,i)+abs(a(j,i))**2
         else if (19<=j.and.j<=20) then
           cu1(k+nt*2,i)=cu1(k+nt*2,i)+abs(a(j,i))**2
         else if (21<=j.and.j<=22) then
           cl1(k+nt*2,i)=cl1(k+nt*2,i)+abs(a(j,i))**2
         else if (23<=j.and.j<=24) then
           cl2(k+nt*2,i)=cl2(k+nt*2,i)+abs(a(j,i))**2
         else if (25<=j.and.j<=28) then
           c12(k+nt*2,i)=c12(k+nt*2,i)+abs(a(j,i))**2
         else
           cu2(k+nt*2,i)=cu2(k+nt*2,i)+abs(a(j,i))**2
         end if
      end do
   end do
end do

a=0
kx=kmax*0.75
ky=kmax*0.75
kz=0.0d0

do k=1,nt
   h=0.75*kmax*(k-1)/nt
   kx=0.75*kmax-h
   ky=0.75*kmax-h

call band_30(kx,ky,kz,a)

call zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)

  do i=1,n
     wc(k+nt*3,i)=w(i)*13.6d0
  end do

   do i=1,n
      do j=1,n
         if (1<=j.and.j<=6) then
           cl25(k+nt*3,i)=cl25(k+nt*3,i)+abs(a(j,i))**2
         else if (7<=j.and.j<=12) then
           c15(k+nt*3,i)=c15(k+nt*3,i)+abs(a(j,i))**2
         else if (13<=j.and.j<=18) then
           cu25(k+nt*3,i)=cu25(k+nt*3,i)+abs(a(j,i))**2
         else if (19<=j.and.j<=20) then
           cu1(k+nt*3,i)=cu1(k+nt*3,i)+abs(a(j,i))**2
         else if (21<=j.and.j<=22) then
           cl1(k+nt*3,i)=cl1(k+nt*3,i)+abs(a(j,i))**2
         else if (23<=j.and.j<=24) then
           cl2(k+nt*3,i)=cl2(k+nt*3,i)+abs(a(j,i))**2
         else if (25<=j.and.j<=28) then
           c12(k+nt*3,i)=c12(k+nt*3,i)+abs(a(j,i))**2
         else
           cu2(k+nt*3,i)=cu2(k+nt*3,i)+abs(a(j,i))**2
         end if
      end do
   end do
end do

end subroutine fullband
