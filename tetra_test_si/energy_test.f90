program test
implicit none

integer::i,k
real(8)::nkx,nky,nkz
real(8)::kx,ky,kz,a,m0,a0,e0,hbar
real(8),dimension(2)::eigen_energy
nky=0
nkz=0

do k=1,51
   nkx=k
call Hami_slove(nkx,nky,nkz,eigen_energy)

do i=1,2
write(1,*) nkx,eigen_energy(i)*13.6
end do

end do

end program test

