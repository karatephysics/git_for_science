program main

implicit none

real(8)::sa,laser_ev
real(8)::kx,ky,kz
integer::i,j,ss,nkx,nky,nkz,atomnum
integer,parameter::n=30
integer,dimension(4)::nt
real(8),dimension(10)::para
real(8),dimension(1000)::energy
real(8),dimension(n,n)::gap
real(8),dimension(n,n)::eigenv
real(8),dimension(n,n)::pump_p
real(8),dimension(10)::dipolepara
real(8)::dipolemoment

!read(*,*) atomnum
atomnum=14

para=0
para(1)=101
nkz=1
gap=0

!read(*,*) laser_ev
laser_ev=2

do i=35,150
kx=i
ky=i+1
energy(i)=i*i
         write(i,"(3f10.4)")  kx,ky,energy(i)
end do

end program main

