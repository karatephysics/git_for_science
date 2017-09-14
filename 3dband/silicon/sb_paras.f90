subroutine paras(atomnum,lattice,energypara,dipolepara)
implicit none 

integer,intent(in)::atomnum
real(8)::lattice
real(8),dimension(20)::energypara,dipolepara

energypara=0
dipolepara=0
!energypara(1)  energypara(1)
!El25 energypara(2) 
!El2 energypara(3)  
!E15 energypara(4)  
!Eu1 energypara(5)  
!E12 energypara(6)  
!Eu25 energypara(7)  
!Eu2 energypara(8)
!so25 energypara(11)
!so15 energypara(12)
!という対称性でエネルギーの低い順、fishmanの論文参照
!10番台はspin-orbit

if (atomnum==14) then
    lattice=5.43d0/0.529177d0  !格子定数/ボーア半径

    energypara(1)=-12.92d0/13.6d0
    energypara(2)=0.0d0
    energypara(3)=4.185d0/13.6d0
    energypara(4)=3.40d0/13.6d0
    energypara(5)=7.07d0/13.6d0
    energypara(6)=9.66d0/13.6d0
    energypara(7)=12.78d0/13.6d0
    energypara(8)=13.46d0/13.6d0
    energypara(11)=0.044d0/13.6d0/3
    !energypara(12)=0.0d0/3

    dipolepara(1)=1.211d0
    dipolepara(2)=0.296d0
    dipolepara(3)=0.542d0
    dipolepara(4)=1.236d0
    dipolepara(5)=1.044d0
    dipolepara(6)=0.742d0
    dipolepara(7)=0.574d0
    dipolepara(8)=0.851d0
    dipolepara(9)=1.097d0
    dipolepara(10)=0.283d0
else if (atomnum==32) then
    lattice=5.6754d0/0.529177d0  !格子定数/ボーア半径

    energypara(1)=-13.14d0/13.6d0
    energypara(2)=0.0d0
    energypara(3)=0.9d0/13.6d0
    energypara(4)=3.22d0/13.6d0
    energypara(5)=7.77d0/13.6d0
    energypara(6)=10.47d0/13.6d0
    energypara(7)=17.0d0/13.6d0
    energypara(8)=18.36d0/13.6d0
    energypara(11)=0.29d0/13.6d0/3
    energypara(12)=0.210d0/13.6d0/3

    dipolepara(1)=1.345d0
    dipolepara(2)=0.019d0
    dipolepara(3)=0.429d0
    dipolepara(4)=1.424d0
    dipolepara(5)=1.139d0
    dipolepara(6)=0.948d0
    dipolepara(7)=0.619d0
    dipolepara(8)=1.076d0
    dipolepara(9)=1.145d0
    dipolepara(10)=0.281d0
else
    lattice=0
    energypara=0
    dipolepara=0
end if 

end subroutine paras

