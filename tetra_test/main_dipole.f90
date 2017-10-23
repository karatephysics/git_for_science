program total_DOS

implicit none
real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
integer::i,j,k,ss,nkx,nky,nkz,indexx
real(8),dimension(4)::DOSpara
real(8),dimension(6)::vol_tetra
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(4,100)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(10001)::DOS,NOS,DOS_parts,di_moment
real(8),dimension(10001)::energy
integer,parameter::N=2000
integer t1, t2, t_rate, t_max, diff


DOSpara(1)=51
!DOSpara(1)=201
!DOSpara(2)=10001
DOSpara(2)=10001
DOSpara(3)=0.8
DOSpara(4)=3

DOS=0
e_ini=DOSpara(3)
e_fin=DOSpara(4)
nn=DOSpara(2)
sa=(e_fin-e_ini)/(nn-1)

!write(*,*) e_ini,sa,100

do ss=1,10001
   energy(ss)=e_ini+sa*(ss-1)
!write(11,*) energy(ss)
end do

do nkz=1,int((int(DOSpara(1))-1)/2)
  do nky=nkz,int(0.75*(int(DOSpara(1))-1)-nkz*0.5)
      do nkx=nky,int(1.5*(int(DOSpara(1))-1)-nkz-nky)!.and.nkx<(DOSpara(1)-1)
!         if (nkx>=(DOSpara(1)-1)) then
!            go to 12
!            else

            call dipole(nkx,nky,nkz,di_moment)

!            end if
         end do
12         sa=sa
      end do
end do

do i=1,10000
   write(1,*) energy(i),DOS(i)
end do

do i=1,10000
   write(2,*) energy(i),NOS(i)
end do

do i=1,10000
   write(3,*) energy(i),di_moment(i)
end do





end program total_DOS

