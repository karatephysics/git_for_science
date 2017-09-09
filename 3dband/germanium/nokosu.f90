program bandgermanium

implicit none

real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
real(8)::kx,ky,kz
real(8)::moji
integer::i,j,ss,nkx,nky,nkz
integer,parameter::n=30
integer,dimension(4)::nt
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(10)::para
real(8),dimension(n)::energy
real(8),dimension(n,n)::gap
real(8),dimension(n,8)::eigenv
real(8)::P,P1,P2,P3,Q,Q1,T,T1,R,R1

para=0
para(1)=101
nkz=1
P=1.345d0
Q=1.139d0
R=0.619d0
P2=0.429d0
P1=0.019d0
Q1=0.948d0
R1=1.076d0
P3=1.424d0
T=1.145d0
T1=0.281d0

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
                i=9
                j=7
                gap(i,j)=energy(i)-energy(j)
                gap(10,10)=0.9d0-energy(j)
                write(18,*) kx,ky,gap(i,j)
                write(19,*) kx,ky,gap(10,10)
                moji=abs(gap(i,j)-1.165d0)
                if (0<=moji.and.moji<=0.01d0) then
                    if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,6)) then
                    !jからiへの遷移を過程
                    para(4)=0
                    para(4)=para(4)+P**2*abs(eigenv(i,1)*eigenv(j,6))
                    para(4)=para(4)+P1**2*abs(eigenv(i,3)*eigenv(j,6))  
                    para(4)=para(4)+P2**2*abs(eigenv(i,1)*eigenv(j,8))
                    para(4)=para(4)+P3**2*abs(eigenv(i,3)*eigenv(j,8))
                    para(4)=para(4)+Q**2*abs(eigenv(i,1)*eigenv(j,2))
                    para(4)=para(4)+Q1**2*abs(eigenv(i,3)*eigenv(j,2))
                    para(4)=para(4)+R**2*abs(eigenv(i,1)*eigenv(j,7))
                    para(4)=para(4)+R1**2*abs(eigenv(i,3)*eigenv(j,7))
                    para(4)=para(4)+T**2*abs(eigenv(i,5)*eigenv(j,2))
                    para(4)=para(4)+T1**2*abs(eigenv(i,4)*eigenv(j,2))
                    write(20,*) kx,ky,gap(i,j),para(4)
                else
                 !   write(20,*) kx,ky,"nan"
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

