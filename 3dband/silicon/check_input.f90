program germanium
implicit none

character, parameter :: jobz="v", uplo="l"
integer, parameter ::n=30, lda=n, lwork=2*n-1
integer ::info, i, j, k
double precision, dimension(3*n-2) :: rwork
double precision, dimension(n) ::w
complex(kind(0d0)) ::a(lda,n), work(lwork)

integer,parameter::nt=1000
real(8)::laser_ev
real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1
real(8),dimension(4*nt,n)::wc
real(8),dimension(4*nt,n)::cl25,cu25,c15,cu1,cl1,cl2,c12,cu2

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

laser_ev=2

call fullband(laser_ev,nt,wc,cl25,c15,cu25,cu1,cl1,cl2,c12,cu2)

do i=1,4*nt
      write(12,*) i-1,wc(i,1)
      write(13,*) i-1,wc(i,3)
      write(15,*) i-1,wc(i,5)
      write(17,*) i-1,wc(i,7)
      write(19,*) i-1,wc(i,9)
      write(111,*) i-1,wc(i,11)
      write(113,*) i-1,wc(i,13)
      write(115,*) i-1,wc(i,15)
      write(117,*) i-1,wc(i,17)
      write(119,*) i-1,wc(i,19)
      write(121,*) i-1,wc(i,21)
      write(123,*) i-1,wc(i,23)
      write(125,*) i-1,wc(i,25)
      write(127,*) i-1,wc(i,27)
      write(129,*) i-1,wc(i,29)
end do

    write(116,"(10A10)") "i-1","cl1","cl25","cl2","c15","cu1","c12","cu25","cu2","total"
do i=1,4*nt
    do j=1,30
        laser_ev=cl25(i,j)+c15(i,j)+cu25(i,j)+cu1(i,j)+cl1(i,j)+cl2(i,j)+c12(i,j)+cu2(i,j)
    end do
        !write(116,"(I10.2,9f10.4)") i-1,cl1(i,30),cl25(i,j),cl2(i,j),c15(i,j),cu1(i,j),c12(i,j),cu25(i,j),cu2(i,j),laser_ev
        !jは30の値になるのでトータルと異なる値になるので注意
        write(116,"(I10.2,9f10.4)") i-1,laser_ev
        laser_ev=0
end do

end program germanium
