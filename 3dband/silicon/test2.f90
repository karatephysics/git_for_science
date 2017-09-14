subroutine bandinfo(kx,ky,laser_ev,energy,dipolepara,eigenv,gap,pump_p,dipolemoment)
implicit none
real(8),intent(in)::kx,ky,laser_ev
real(8),dimension(30),intent(in)::energy
real(8),dimension(10),intent(in)::dipolepara
real(8),dimension(30,30),intent(in)::eigenv
integer::i,j
real(8)::dipolemoment
real(8),dimension(10)::para
real(8),dimension(30,30)::pump_p
real(8),dimension(30,30)::gap

para=0
gap=0
pump_p=0

!全バンド計算がしたい場合はdoループで行う
!do i=9,30
!    do j=1,8
i=9
j=8
!gap(i,j)=energy(i)-energy(j)
gap(10,10)=1.12-energy(j)
!write(18,*) kx,ky,gap(i,j)
write(18,*) kx,ky,gap(10,10)
!para(4)=abs(gap(i,j)-laser_ev)
para(1)=abs(gap(10,10)-laser_ev)
if (0<=para(1).and.para(1)<=0.01) then
    !ギャップは+-6-7nmとする,para(4)<=0.01の場合、para(1)を小さくしないとサーチにひっかからない
    pump_p(i,j)=1
    para(2)=0
    !jからiへの遷移を過程
    if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,6)) then
    para(2)=para(2)+dipolepara(1)*abs(eigenv(i,1)*eigenv(j,6)+eigenv(j,1)*eigenv(i,6))
    if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,6)) then
    para(2)=para(2)+dipolepara(2)*abs(eigenv(i,3)*eigenv(j,6)+eigenv(j,3)*eigenv(i,6)) 
    if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,8)) then
    para(2)=para(2)+dipolepara(3)*abs(eigenv(i,1)*eigenv(j,8)+eigenv(j,1)*eigenv(i,8))  
    if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,8)) then
    para(2)=para(2)+dipolepara(4)*abs(eigenv(i,3)*eigenv(j,8)+eigenv(j,3)*eigenv(i,8))  
    if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,2)) then
    para(2)=para(2)+dipolepara(5)*abs(eigenv(i,1)*eigenv(j,2)+eigenv(j,1)*eigenv(i,2)) 
    if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,2)) then
    para(2)=para(2)+dipolepara(6)*abs(eigenv(i,3)*eigenv(j,2)+eigenv(j,3)*eigenv(i,2))  
    if (0.000<=eigenv(i,1).and.0.000<=eigenv(j,7)) then
    para(2)=para(2)+dipolepara(7)*abs(eigenv(i,1)*eigenv(j,7)+eigenv(j,1)*eigenv(i,7)) 
    if (0.000<=eigenv(i,3).and.0.000<=eigenv(j,7)) then
    para(2)=para(2)+dipolepara(8)*abs(eigenv(i,3)*eigenv(j,7)+eigenv(j,3)*eigenv(i,7)) 
    if (0.000<=eigenv(i,5).and.0.000<=eigenv(j,2)) then
    para(2)=para(2)+dipolepara(9)*abs(eigenv(i,5)*eigenv(j,2)+eigenv(j,5)*eigenv(i,2)) 
    if (0.000<=eigenv(i,4).and.0.000<=eigenv(j,2)) then
    para(2)=para(2)+dipolepara(10)*abs(eigenv(i,4)*eigenv(j,2)+eigenv(j,4)*eigenv(i,2)) 
    dipolemoment=para(2)**2
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
    !双極子遷移モーメントのあたいはpara(2)
    write(211,"(2f10.4,2i10,f10.4)") gap(10,10),para(1),i,j,dipolemoment
    write(20,*) kx,ky,gap(10,10)
else
    write(20,*) kx,ky,"nan"
end if
!end do
!end do

end subroutine bandinfo

