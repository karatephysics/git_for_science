program main

implicit none

real(8)::sa,laser_ev
real(8)::kx,ky,kz
integer::i,j,ss,nkx,nky,nkz,atomnum
integer,parameter::n=30
integer,dimension(4)::nt
real(8),dimension(10)::para
real(8),dimension(n)::energy
real(8),dimension(n,n)::gap
real(8),dimension(n,n)::eigenv
real(8),dimension(n,n)::pump_p
real(8),dimension(10)::dipolepara
real(8)::dipolemoment

read(*,*) atomnum

para=0
para(1)=101
nkz=1
gap=0

!read(*,*) laser_ev
laser_ev=2

write(211,"(5a10)") 'laser_ev','diff','i','j','dipole'

do nky=nkz,int(para(1))
    do nkx=nkz,int(para(1))
        para(2)=int(0.75*int((para(1)-1))-nkz*0.5+1)
        if (nky<=int(para(2))) then
            para(3)=int(1.5*int((para(1)-1))-nkz-nky+1)
            if (nky<=nkx.and.nkx<=int((para(1)-1)).and.nkx<=para(3)) then
                call ham(atomnum,nkx-1,nky-1,nkz-1,int(para(1)),energy,eigenv)
                kx=real((nkx-1))/para(1)!100
                ky=real((nky-1))/para(1)!100

                write(7,*)  kx,ky,energy(3)
                write(9,*)  kx,ky,energy(5)
                write(11,*)  kx,ky,energy(7)
                write(13,*)  kx,ky,energy(9)
                write(15,*)  kx,ky,energy(11)
                write(17,*)  kx,ky,energy(13)
                call bandinfo(kx,ky,laser_ev,energy,dipolepara,eigenv,gap,pump_p,dipolemoment)
!                write(211,"(2f10.4)") gap(10,10),dipolemoment
            else
                kx=real((nkx-1))/para(1)!100
                ky=real((nky-1))/para(1)!100
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
            kx=real((nkx-1))/para(1)!100
            ky=real((nky-1))/para(1)!100
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

end program main

