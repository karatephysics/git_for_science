program bandge

implicit none
real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
real(8)::kx,ky,kz
real(8)::moji
integer::i,j,ss,nkx,nky,nkz,para1,para2
real(8),dimension(4)::DOSpara
real(8),dimension(6)::vol_tetra
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(4,2)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(10001)::DOS,NOS,DOS_parts
real(8),dimension(30)::energy
real(8),dimension(30,30)::gap

DOSpara(1)=101

nkz=1

do nky=nkz,int(DOSpara(1))
    do nkx=nkz,int(DOSpara(1))
        para1=int(0.75*(int(DOSpara(1))-1)-nkz*0.5+1)
        if (nky<=para1) then
        para1=int(1.5*(int(DOSpara(1))-1)-nkz-nky+1)
        if (nky<=nkx.and.nkx<=(DOSpara(1)-1).and.nkx<=para1) then

        call ham(nkx-1,nky-1,nkz-1,energy)
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
                gap(10,10)=1.12d0-energy(j)
                write(18,*) kx,ky,gap(i,j)
                write(19,*) kx,ky,gap(10,10)
                moji=abs(gap(10,10)-1.55d0)
                if (0<=moji.and.moji<=0.01d0) then
                  write(20,*) kx,ky,gap(10,10)
                else
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

end program bandge

