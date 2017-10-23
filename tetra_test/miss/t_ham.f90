program test
implicit none
integer::i
real(8)::a,b,c
real(8),dimension(2)::energy

a=1
b=0
c=0
call ham(a,b,c,energy)

do i=1,2
   write(*,'(f10.8)') energy(i)*13.6d0
end do

end program test


subroutine ham(a,b,c,energy)
implicit none
real(8),intent(in)::a,b,c
real(8),dimension(2)::energy
real(8),dimension(2)::DOSpara
real(8)::kx,ky,kz,a0,egap,pi

DOSpara(1)=51

energy=0
egap=0.5d0
pi=3.141592d0
a0=5.43d0/0.529d0

kx=(pi/a0)*a/(DOSpara(1)-1)
ky=(pi/a0)*b/(DOSpara(1)-1)
kz=(pi/a0)*c/(DOSpara(1)-1)

energy(1)=-(kx**2+ky**2+kz**2)
energy(2)=kx**2+ky**2+kz**2+egap/13.6

end subroutine ham

