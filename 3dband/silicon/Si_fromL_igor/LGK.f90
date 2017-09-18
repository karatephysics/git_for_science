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
do nkz=1,int(para(1))
    do nky=1,int(para(1))
                nkx=nky
        para(2)=int(3*para(1)/4+1-nkz*0.5)
        if (nky<=int(para(2))) then
            para(3)=int(1.5*para(1)+2-nkz-nky)
            if (nkx<=nky.and.nkx<=int(para(1)).and.nkx<=para(3).and.nkz<=nky) then
                    write(227,"(3I10.3,2f10.3)") nkx,nky,nkz,para(2),para(3)!,nkx+nky+nkz,int(1.5*para(1))
                    !write(*,"(4I10.3)") nkx,nky,nkz,nkx+nky+nkz
!                if (nkx+nky+nkz==int(1.5*para(1)+1)) then  !L-G-X
!               if (nky==nkz) then  !L-K-W
!              if (nkx==nky) then  !L-G-K
!                if (nkx==int(para(1)+1)) then  !U-X-W
!                if (nkz==1) then  !G-X-K
                    write(226,"(5I10.3)") nkx-1,nky-1,nkz-1!,nkx+nky+nkz,int(1.5*para(1))
                call ham(atomnum,nkx-1,nky-1,nkz-1,int(para(1)),energy,eigenv)
                kx=real((nkx-1))/100
                ky=real((nky-1))/100
                kz=real((nkz-1))/100
moji=0.75*kx+0.75*ky
                write(7,*)  kz,moji,energy(3)
                write(9,*)  kz,moji,energy(5)
                write(11,*)  kz,moji,energy(7)
                write(13,*)  kz,moji,energy(9)
                write(15,*)  kz,moji,energy(11)
                write(17,*)  kz,moji,energy(13)
!                call bandinfo(kx,ky,laser_ev,energy,dipolepara,eigenv,gap,pump_p,dipolemoment)
!                end if
                i=9
                j=7
                gap(i,j)=energy(i)-energy(j)
                gap(10,10)=1.12d0-energy(j)
                write(18,*) kz,moji,gap(i,j)
                write(19,*) kz,moji,gap(10,10)
                moji=abs(gap(10,10)-1.55d0)
                if (0<=moji.and.moji<=0.01d0) then
                  write(20,*) kz,moji,gap(10,10)
                else
                  write(20,*) kz,moji,"nan"
                end if
!                write(211,"(2f10.4)") gap(10,10),dipolemoment
            else
                kx=real((nkx-1))/100
                ky=real((nky-1))/100
                kz=real((nkz-1))/100
moji=0.75*kx+0.75*ky
                write(7,*)  kz,moji,"nan"
                write(9,*)  kz,moji,"nan"
                write(11,*)  kz,moji,"nan"
                write(13,*)  kz,moji,"nan"
                write(15,*)  kz,moji,"nan"
                write(17,*)  kz,moji,"nan"
                write(18,*) kz,moji,"nan"
                write(19,*) kz,moji,"nan"
                write(20,*) kz,moji,"nan"
            end if
        else
                kx=real((nkx-1))/100
                ky=real((nky-1))/100
                kz=real((nkz-1))/100
moji=0.75*kx+0.75*ky
                write(7,*)  kz,moji,"nan"
                write(9,*)  kz,moji,"nan"
                write(11,*)  kz,moji,"nan"
                write(13,*)  kz,moji,"nan"
                write(15,*)  kz,moji,"nan"
                write(17,*)  kz,moji,"nan"
                write(18,*) kz,moji,"nan"
                write(19,*) kz,moji,"nan"
                write(20,*) kz,moji,"nan"
        end if
    11       sa=sa
    end do
end do
!end do

end program main

