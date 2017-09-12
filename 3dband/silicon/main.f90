program bandgermanium

implicit none

real(8)::sa,laser_ev
real(8)::kx,ky,kz
integer::i,j,ss,nkx,nky,nkz
integer,parameter::n=30
integer,dimension(4)::nt
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(10)::para
real(8),dimension(n)::energy
real(8),dimension(n,n)::gap
real(8),dimension(n,n)::eigenv
real(8),dimension(n,n)::pump_p
real(8),dimension(10)::dipolepara
real(8)::dipolemoment

para=0
para(1)=201
nkz=1
gap=0
dipolepara(1)=1.211d0 !P
dipolepara(2)=0.296d0 !P1
dipolepara(3)=0.542d0 !P2
dipolepara(4)=1.236d0 !P3
dipolepara(5)=1.044d0 !Q
dipolepara(6)=0.742d0 !Q1
dipolepara(7)=0.574d0 !R
dipolepara(8)=0.851d0 !R1
dipolepara(9)=1.097d0 !T
dipolepara(10)=0.283d0 !T1

read(*,*) laser_ev

write(211,"(5a10)") 'laser_ev','diff','i','j','dipole'

do nky=nkz,int(para(1))
    do nkx=nkz,int(para(1))
        para(2)=int(0.75*int((para(1)-1))-nkz*0.5+1)
        if (nky<=int(para(2))) then
            para(3)=int(1.5*int((para(1)-1))-nkz-nky+1)
            if (nky<=nkx.and.nkx<=int((para(1)-1)).and.nkx<=para(3)) then
                call ham(nkx-1,nky-1,nkz-1,int(para(1)),energy,eigenv)
                kx=real((nkx-1))/100
                ky=real((nky-1))/100

                write(7,*)  kx,ky,energy(3)
                write(9,*)  kx,ky,energy(5)
                write(11,*)  kx,ky,energy(7)
                write(13,*)  kx,ky,energy(9)
                write(15,*)  kx,ky,energy(11)
                write(17,*)  kx,ky,energy(13)
                call bandinfo(kx,ky,laser_ev,energy,dipolepara,eigenv,gap,pump_p,dipolemoment)
!                write(211,"(2f10.4)") gap(10,10),dipolemoment
            else
                kx=real((nkx-1))/100
                ky=real((nky-1))/100
                write(7,*)  kx,ky,"nan"
                write(9,*)  kx,ky,"nan"
                write(11,*)  kx,ky,"nan"
                write(13,*)  kx,ky,"nan"
                write(15,*)  kx,ky,"nan"
                write(17,*)  kx,ky,"nan"
                write(18,*) kx,ky,"nan"
                write(19,*) kx,ky,"nan"
                write(20,*) kx,ky,"nan"
            end if
        else
            kx=real((nkx-1))/100
            ky=real((nky-1))/100
            write(7,*) kx,ky,"nan"
            write(9,*) kx,ky,"nan"
            write(11,*) kx,ky,"nan"
            write(13,*) kx,ky,"nan"
            write(15,*) kx,ky,"nan"
            write(17,*) kx,ky,"nan"
            write(18,*) kx,ky,"nan"
            write(19,*) kx,ky,"nan"
            write(20,*) kx,ky,"nan"
        end if
    11       sa=sa
    end do
end do

end program bandgermanium

