program total_DOS

implicit none
real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
integer::i,j,ss,nkx,nky,nkz,indexx
real(8),dimension(4)::DOSpara
real(8),dimension(6)::vol_tetra
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(4,2)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(10001)::DOS,NOS,DOS_parts
real(8),dimension(10001)::energy

DOSpara(1)=51
DOSpara(2)=10001
DOSpara(3)=-1
DOSpara(4)=2

DOS=0
e_ini=DOSpara(3)
e_fin=DOSpara(4)
nn=DOSpara(2)
sa=(e_fin-e_ini)/(nn-1)

do ss=1,10001
   energy(ss)=e_ini+sa*(ss-1)
end do

nkz=2
nky=4
nkx=18

do nkz=1,int((int(DOSpara(1))-1)*1/2)
   do nky=nkz,int(3/4*(int(DOSpara(1))-1)-nkz/2)
         write(10,*) nky
      do nkx=nky,int(3/2*(int(DOSpara(1))-1)-nkz-nky)!.and.nkx<(DOSpara(1)-1)

         call tetra_apex_output(nkx,nky,nkz,tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,vol_tetra)
         do i=1,3
            do j=1,4
               write(26,*) i,j,tetra1(i,j)
            end do
         end do

         call tetra_energy_output(tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,&
         ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6)

         do i=1,4
            do j=1,2
               write(27,*) i,j,ene_tetra1(i,j)
            end do
         end do

            do indexx=1,2
               call parts_cal(ene_tetra1(1,indexx),ene_tetra1(2,indexx),&
               ene_tetra1(3,indexx),ene_tetra1(4,indexx),vol_tetra(1),DOS)
               call parts_cal(ene_tetra2(1,indexx), ene_tetra2(2,indexx),&
               ene_tetra2(3,indexx),ene_tetra2(4,indexx),vol_tetra(2),DOS)
               call parts_cal(ene_tetra3(1,indexx), ene_tetra3(2,indexx),&
               ene_tetra3(3,indexx),ene_tetra3(4,indexx),vol_tetra(3),DOS)
               call parts_cal(ene_tetra4(1,indexx), ene_tetra4(2,indexx),&
               ene_tetra4(3,indexx),ene_tetra4(4,indexx),vol_tetra(4),DOS)
               call parts_cal(ene_tetra5(1,indexx), ene_tetra5(2,indexx),&
               ene_tetra5(3,indexx),ene_tetra5(4,indexx),vol_tetra(5),DOS)
               call parts_cal(ene_tetra6(1,indexx), ene_tetra6(2,indexx),&
               ene_tetra6(3,indexx),ene_tetra6(4,indexx),vol_tetra(6),DOS)
              ! stock=stock+1
              do i=1,10000
                 write(28,*) DOS(i)
              end do

            end do
         end do
      end do
end do

do i=1,10000
   write(1,*) energy(i),DOS(i)
end do

end program total_DOS
