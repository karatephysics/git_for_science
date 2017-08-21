program germanium
implicit none

character, parameter :: jobz="v", uplo="l"
integer, parameter ::n=30, lda=n, lwork=2*n-1
integer ::info, i, j, k
double precision, dimension(3*n-2) :: rwork
double precision, dimension(n) ::w
complex(kind(0d0)) ::a(lda,n), work(lwork)

integer,parameter::nt=100
real(8)::laser_ev,den
real(8)::x,l,g,g2,g3,s
real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1
real(8),dimension(4*nt,n,n)::gap !エネルギーギャップ
real(8),dimension(4*nt,n)::wc
real(8),dimension(4*nt,n)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2
real(8),dimension(4*nt,n)::deri
real(8),dimension(4*nt,n,n)::density
real(8),dimension(4*nt,n,n)::pump_p,e

!P=1.345d0
!Q=1.139d0
!R=0.619d0
!P2=0.429d0
!P1=0.019d0
!Q1=0.948d0
!R1=1.076d0
!P3=1.424d0
!T=1.145d0
!T1=0.281d0

read(*,*) laser_ev

call fullband(laser_ev,nt,wc,gap,density,pump_p,cl25,c15,cu25,cu1,cl1,cl2,c12,cu2)

do i=1,4*nt
      write(12,*) i-1,wc(i,1)
      write(13,*) i-1,wc(i,3)
      write(15,*) i-1,wc(i,5)
      write(17,*) i-1,wc(i,7)
      write(19,*) i-1,wc(i,9)
      write(111,*) i-1,wc(i,11)
      write(113,*) i-1,wc(i,13)
      write(115,*) i-1,wc(i,15)
end do


end program germanium
