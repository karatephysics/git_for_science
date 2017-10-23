subroutine Hami_slove(nkx,nky,nkz,eigen_energy)

real(8),intent(in)::nkx,nky,nkz
real(8)::kx,ky,kz,a,m0,a0,e0,hbar,egap
real(8),dimension(4)::DOSpara
real(8),dimension(2)::eigen_energy

DOSpara(1)=51
!m0=9.109e-31//kg
!a0=5.43e-10//m
!e0=1.602e-19//C
!hbar=1.054e-34//J*sec
m0=1
pi=3.141592
a0=5.658/0.529
hbar=1

mxx=m0*1
myy=m0*1
mzz=m0*1
myz=m0*1
mzx=m0*1
mxy=m0*1
kx=(2*pi/a0)*nkx/(DOSpara(1)-1)
ky=(2*pi/a0)*nky/(DOSpara(1)-1)
kz=(2*pi/a0)*nkz/(DOSpara(1)-1)
eigen_energy=0
egap=0.5

eigen_energy(2)=kx**2+ky**2+kz**2+egap/13.6
!eigen_energy(2)=(hbar**2/(2*mxx))*kx**2+(hbar**2/(2*myy))*ky**2+(hbar**2/(2*mzz))*kz**2
!eigen_energy(2)=eigen_energy(2)*2
!eigen_energy(2)=eigen_energy(2)+e0/2
!eigen_energy(1)=eigen_energy(1)-(kx**2+ky**2+*kz**2)
eigen_energy(1)=-(kx**2+ky**2+kz**2)
!eigen_energy(1)=eigen_energy(1)*1
!eigen_energy=eigen_energy/e0

end subroutine Hami_slove
