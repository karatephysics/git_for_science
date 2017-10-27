program total_DOS

implicit none
real(8)::e_ini,e_fin,nn,sa,e0,e1,e2,e3
integer::i,j,k,ss,nkx,nky,nkz,indexx
real(8),dimension(4)::para
real(8),dimension(6)::vol_tetra
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(4,100)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(10001)::DOS,NOS,DOS_parts
real(8),dimension(10001)::energy
!integer,parameter::N=2000
integer t1, t2, t_rate, t_max, diff


para(1)=201
!para(1)=201
!para(2)=10001
para(2)=2001
para(3)=0.8d0
para(4)=10.0d0

DOS=0
NOS=0
e_ini=para(3)
e_fin=para(4)
nn=para(2)
sa=(e_fin-e_ini)/(nn-1)

!write(*,*) e_ini,sa,100

do ss=1,int(para(2))
   energy(ss)=e_ini+sa*(ss-1)
!write(11,*) energy(ss)
end do

do nkz=1,int((int(para(1))-1)/2)
call system_clock(t1)   ! 開始時を記録
  do nky=nkz,int(0.75*(int(para(1))-1)-nkz*0.5)
      do nkx=nky,int(1.5*(int(para(1))-1)-nkz-nky)!.and.nkx<(para(1)-1)
         if (nkx>=(para(1)-1)) then
            go to 11
            else

         call tetra_apex_output(para,nkx,nky,nkz,tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,vol_tetra)

         call tetra_energy_output(tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,&
         ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6)

            do indexx=1,44
               call parts_cal(para,ene_tetra1(1,indexx),ene_tetra1(2,indexx),&
               ene_tetra1(3,indexx),ene_tetra1(4,indexx),vol_tetra(1),DOS,NOS)
               call parts_cal(para,ene_tetra2(1,indexx), ene_tetra2(2,indexx),&
               ene_tetra2(3,indexx),ene_tetra2(4,indexx),vol_tetra(2),DOS,NOS)
               call parts_cal(para,ene_tetra3(1,indexx), ene_tetra3(2,indexx),&
               ene_tetra3(3,indexx),ene_tetra3(4,indexx),vol_tetra(3),DOS,NOS)
               call parts_cal(para,ene_tetra4(1,indexx), ene_tetra4(2,indexx),&
               ene_tetra4(3,indexx),ene_tetra4(4,indexx),vol_tetra(4),DOS,NOS)
               call parts_cal(para,ene_tetra5(1,indexx), ene_tetra5(2,indexx),&
               ene_tetra5(3,indexx),ene_tetra5(4,indexx),vol_tetra(5),DOS,NOS)
               call parts_cal(para,ene_tetra6(1,indexx), ene_tetra6(2,indexx),&
               ene_tetra6(3,indexx),ene_tetra6(4,indexx),vol_tetra(6),DOS,NOS)
              ! stock=stock+1
            end do
            end if
         end do
11         sa=sa
      end do
  call system_clock(t2, t_rate, t_max)   ! 終了時を記録
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  print "(I5, A, F10.3)", nkz,"time it took was:",  diff/dble(t_rate)

end do

do i=1,int(para(2))
   write(1,*) energy(i),DOS(i)
end do

do i=1,int(para(2))
   NOS(i)=NOS(i)/10000000
end do

do i=1,int(para(2))
   write(2,*) energy(i),NOS(i)
end do



end program total_DOS

