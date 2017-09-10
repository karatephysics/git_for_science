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
real(8),dimension(n,8)::eigenv
real(8),dimension(n,8)::pump_p
real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1

para=0
para(1)=201
nkz=1
gap=0
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

!                全部のバンドの計算がしたい場合はdoループを解除
  !              do i=9,30
  !                  do j=1,8
                        i=9
                        j=8
                        !gap(i,j)=energy(i)-energy(j)
                        gap(10,10)=1.12-energy(j)
                        !write(18,*) kx,ky,gap(i,j)
                        write(18,*) kx,ky,gap(10,10)
                        !para(4)=abs(gap(i,j)-laser_ev)
                        para(4)=abs(gap(10,10)-laser_ev)
                        if (0<=para(4).and.para(4)<=0.01) then
                            !ギャップは+-6-7nmとする,para(4)<=0.01の場合、para(1)を小さくしないとサーチにひっかからない
                            pump_p(i,j)=1
                            para(5)=0
                            !jからiへの遷移を過程
                            if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,6)) then
                            para(5)=para(5)+P*abs(eigenv(i,1)*eigenv(j,6)+eigenv(j,1)*eigenv(i,6))
                            if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,6)) then
                            para(5)=para(5)+P1*abs(eigenv(i,3)*eigenv(j,6)+eigenv(j,3)*eigenv(i,6)) 
                            if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,8)) then
                            para(5)=para(5)+P2*abs(eigenv(i,1)*eigenv(j,8)+eigenv(j,1)*eigenv(i,8))  
                            if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,8)) then
                            para(5)=para(5)+P3*abs(eigenv(i,3)*eigenv(j,8)+eigenv(j,3)*eigenv(i,8))  
                            if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,2)) then
                            para(5)=para(5)+Q*abs(eigenv(i,1)*eigenv(j,2)+eigenv(j,1)*eigenv(i,2)) 
                            if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,2)) then
                            para(5)=para(5)+Q1*abs(eigenv(i,3)*eigenv(j,2)+eigenv(j,3)*eigenv(i,2))  
                            if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,7)) then
                            para(5)=para(5)+R*abs(eigenv(i,1)*eigenv(j,7)+eigenv(j,1)*eigenv(i,7)) 
                            if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,7)) then
                            para(5)=para(5)+R1*abs(eigenv(i,3)*eigenv(j,7)+eigenv(j,3)*eigenv(i,7)) 
                            if (0.000<=eigenv(i,5).and.0.000<=eigenv(j,2)) then
                            para(5)=para(5)+T*abs(eigenv(i,5)*eigenv(j,2)+eigenv(j,5)*eigenv(i,2)) 
                            if (0.000<=eigenv(i,4).and.0.000<=eigenv(j,2)) then
                            para(5)=para(5)+T1*abs(eigenv(i,4)*eigenv(j,2)+eigenv(j,4)*eigenv(i,2)) 
                            para(5)=para(5)**2
                            end if
                            end if
                            end if
                            end if
                            end if
                            end if
                            end if
                            end if
                            end if
                            end if
                            !双極子遷移モーメントのあたいはpara(5)
                            write(211,"(2f10.4,2i10,f10.4)") gap(10,10),para(4),i,j,para(5)
                            !write(20,*) kx,ky,gap(i,j)
                            write(20,*) kx,ky,gap(10,10)
                        else
                            write(20,*) kx,ky,"nan"
                        end if
    !                end do
   !             end do
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
    11       sa=sa
    end do
end do

end program bandgermanium

