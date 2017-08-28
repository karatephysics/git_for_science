program bandge

implicit none
real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
integer::i,j,ss,nkx,nky,nkz,indexx
real(8),dimension(4)::DOSpara
real(8),dimension(6)::vol_tetra
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(4,2)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(10001)::DOS,NOS,DOS_parts
real(8),dimension(2)::energy
real(8),dimension(30,30)::gap

DOSpara(1)=101

!do nkz=1,int((int(DOSpara(1))-1)/2)
nkz=1
  do nky=nkz,int(0.75*(int(DOSpara(1))-1)-nkz*0.5)
      do nkx=nky,int(1.5*(int(DOSpara(1))-1)-nkz-nky)!.and.nkx<(DOSpara(1)-1)
         if (nkx>=(DOSpara(1)-1)) then
            go to 11
            else
                    call ham(nkx-1,nky-1,nkz-1,energy)
                write(15,*) nkx-1,nky-1,energy(1)*13.6d0
                !write(15,*) nkx-1,nky-1,energy(5)
                write(7,*) nkx-1,nky-1,energy(2)*13.6d0
                !write(7,*) nkx-1,nky-1,energy(7)
                !do i=9,18
                !    do j=1,8
                i=1
                j=2
                    gap(i,j)=energy(i)-energy(j)
                write(17,*) nkx,nky,gap(i,j)
                !    moji=abs(gap(k,i,j)-laser)
                !    if (0<=moji.and.moji<=0.001d0) then
                !       pump_p(k,i,j)=1
                !    end if
                !    gap(k,i,j)=(energy(i)-energy(j))*13.6d0
                !    end do
                !end do
         end if
         end do
11       sa=sa
      end do
!end do

end program bandge

