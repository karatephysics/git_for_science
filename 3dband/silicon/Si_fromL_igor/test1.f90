program main

implicit none

real(8)::sa,laser_ev
real(8)::kx,ky,kz,moji
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

!read(*,*) atomnum
atomnum=14

para=0
para(1)=101
nkz=1
gap=0

!read(*,*) laser_ev
laser_ev=2

write(211,"(5a10)") 'laser_ev','diff','i','j','dipole'

!do nkz=1,int(para(1)+1)
do nky=1,int(para(1))
    do nkx=1,int(para(1))
                nkz=nky
        para(2)=int(3*para(1)/4-nkx*0.5+1)
        if (nky<=int(para(2))) then
            para(3)=int(1.5*para(1)-nkz-nky+2)
            if (nky<=nkx.and.nkx<=int(para(1)).and.nkx<=int(para(3))) then
                write(227,"(3I10.4,2f10.4)") nkx,nky,nkz,para(2),para(3)
                    write(226,"(5I10.3)") nkx-1,nky-1,nkz-1!,nkx+nky+nkz,int(1.5*para(1))
                call ham(atomnum,nkx-1,nky-1,nkz-1,int(para(1)),energy,eigenv)
                kx=real((nkx-1))/para(1)!100
                ky=real((nky-1))/para(1)!100

                do i=101,115
                j=2*(i-100)
                write(i,"(3f10.4)")  kx,ky,energy(j)
                end do
                i=9
                j=7
                gap(i,j)=energy(i)-energy(j)
                gap(10,10)=1.12d0-energy(j)
                write(116,*) kx,ky,gap(i,j)
                write(117,*) kx,ky,gap(10,10)
                moji=abs(gap(10,10)-1.55d0)
                if (0<=moji.and.moji<=0.01d0) then
                  write(118,*) kx,ky,gap(10,10)
                else
                  write(118,*) kx,ky,"nan"
                end if
            else
                kx=real((nkx-1))/para(1)!100
                ky=real((nky-1))/para(1)!100
                do i=101,115
                j=2*(i-100)
                write(i,"(2f10.4,a10)")  kx,ky,"nan"
                end do
                write(116,*) kx,ky,"nan"
                write(117,*) kx,ky,"nan"
                write(118,*) kx,ky,"nan"
            end if
        else
            kx=real((nkx-1))/para(1)!100
            ky=real((nky-1))/para(1)!100
                do i=101,115
                j=2*(i-100)
                write(i,"(2f10.4,a10)")  kx,ky,"nan"
                end do
                write(116,*) kx,ky,"nan"
                write(117,*) kx,ky,"nan"
                write(118,*) kx,ky,"nan"
        end if
    11       sa=sa
    end do
end do
!end do

end program main

