subroutine tetra_apex_output(para,nkx,nky,nkz,tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,vol_tetra)

integer,intent(in)::nkx,nky,nkz
real(8)::nnkx,nnky,nnkz
integer::jj,ss,i
real(8)::nn,det
real(8),dimension(4),intent(in)::para
real(8),dimension(3)::wave0
real(8),dimension(6)::vol_tetra
real(8),dimension(4,8)::cube
real(8),dimension(3,4)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
real(8),dimension(3,3)::mat

nn=para(1)

nnkx=nkx/(nn-1)
nnky=nky/(nn-1)
nnkz=nkz/(nn-1)

wave0(1)=nnkx
wave0(2)=nnky
wave0(3)=nnkz

!do i=1,3
!   write(10,'(3f8.4)') wave0(i)
!end do

tetra1=0
tetra2=0
tetra3=0
tetra4=0
tetra5=0
tetra6=0

if (wave0(3)>1.and.(wave0(1)+wave0(2)+wave0(3))>3/2) then 
   write(*,*) "input point is not in IBZ"
   cube=0
   tetra1=0
   tetra2=0
   tetra3=0
   tetra4=0
   tetra5=0
   tetra6=0
end if

cube=0
cube(1,1)=nnkx
cube(2,1)=nnky
cube(3,1)=nnkz
cube(1,2)=nnkx+1/(nn-1)
cube(2,2)=nnky
cube(3,2)=nnkz
cube(1,3)=nnkx
cube(2,3)=nnky+1/(nn-1)
cube(3,3)=nnkz
cube(1,4)=nnkx+1/(nn-1)
cube(2,4)=nnky+1/(nn-1)
cube(3,4)=nnkz
cube(1,5)=nnkx
cube(2,5)=nnky
cube(3,5)=nnkz+1/(nn-1)
cube(1,6)=nnkx+1/(nn-1)
cube(2,6)=nnky
cube(3,6)=nnkz+1/(nn-1)
cube(1,7)=nnkx
cube(2,7)=nnky+1/(nn-1)
cube(3,7)=nnkz+1/(nn-1)
cube(1,8)=nnkx+1/(nn-1)
cube(2,8)=nnky+1/(nn-1)
cube(3,8)=nnkz+1/(nn-1)

!do i=1,3
!   do j=1,8
!      write(10,'(I6,I6,f8.4)') i,j,cube(i,j)
!   end do
!end do

do jj=1,8
   if (cube(1,jj)>1) then
      cube(4,jj)=1
   end if 
   if (cube(3,jj)>cube(1,jj)) then
      cube(4,jj)=1
   end if 
   if (cube(2,jj)>cube(1,jj)) then
      cube(4,jj)=1
   end if 
   if (cube(3,jj)>cube(2,jj)) then
      cube(4,jj)=1
   end if 
   if ((cube(1,jj)+cube(2,jj)+cube(3,jj))>3/2) then
   cube(3,jj)=1
   end if
end do

tetra1(1,1)=nnkx*(nn-1)
tetra1(2,1)=nnky*(nn-1)
tetra1(3,1)=nnkz*(nn-1)
tetra1(1,2)=nnkx*(nn-1)+1
tetra1(2,2)=nnky*(nn-1)
tetra1(3,2)=nnkz*(nn-1)
tetra1(1,3)=nnkx*(nn-1)+1
tetra1(2,3)=nnky*(nn-1)+1
tetra1(3,3)=nnkz*(nn-1)
tetra1(1,4)=nnkx*(nn-1)+1
tetra1(2,4)=nnky*(nn-1)+1
tetra1(3,4)=nnkz*(nn-1)+1
            
tetra2(1,1)=nnkx*(nn-1)
tetra2(2,1)=nnky*(nn-1)
tetra2(3,1)=nnkz*(nn-1)
tetra2(1,2)=nnkx*(nn-1)+1
tetra2(2,2)=nnky*(nn-1)
tetra2(3,2)=nnkz*(nn-1)
tetra2(1,3)=nnkx*(nn-1)+1
tetra2(2,3)=nnky*(nn-1)
tetra2(3,3)=nnkz*(nn-1)+1
tetra2(1,4)=nnkx*(nn-1)+1
tetra2(2,4)=nnky*(nn-1)+1
tetra2(3,4)=nnkz*(nn-1)+1

tetra3(1,1)=nnkx*(nn-1)
tetra3(2,1)=nnky*(nn-1)
tetra3(3,1)=nnkz*(nn-1)
tetra3(1,2)=nnkx*(nn-1)
tetra3(2,2)=nnky*(nn-1)+1
tetra3(3,2)=nnkz*(nn-1)
tetra3(1,3)=nnkx*(nn-1)+1
tetra3(2,3)=nnky*(nn-1)+1
tetra3(3,3)=nnkz*(nn-1)
tetra3(1,4)=nnkx*(nn-1)+1
tetra3(2,4)=nnky*(nn-1)+1
tetra3(3,4)=nnkz*(nn-1)+1
             
tetra4(1,1)=nnkx*(nn-1)
tetra4(2,1)=nnky*(nn-1)
tetra4(3,1)=nnkz*(nn-1)
tetra4(1,2)=nnkx*(nn-1)
tetra4(2,2)=nnky*(nn-1)+1
tetra4(3,2)=nnkz*(nn-1)
tetra4(1,3)=nnkx*(nn-1)
tetra4(2,3)=nnky*(nn-1)+1
tetra4(3,3)=nnkz*(nn-1)+1
tetra4(1,4)=nnkx*(nn-1)+1
tetra4(2,4)=nnky*(nn-1)+1
tetra4(3,4)=nnkz*(nn-1)+1
             
tetra5(1,1)=nnkx*(nn-1)
tetra5(2,1)=nnky*(nn-1)
tetra5(3,1)=nnkz*(nn-1)
tetra5(1,2)=nnkx*(nn-1)
tetra5(2,2)=nnky*(nn-1)
tetra5(3,2)=nnkz*(nn-1)+1
tetra5(1,3)=nnkx*(nn-1)+1
tetra5(2,3)=nnky*(nn-1)
tetra5(3,3)=nnkz*(nn-1)+1
tetra5(1,4)=nnkx*(nn-1)+1
tetra5(2,4)=nnky*(nn-1)+1
tetra5(3,4)=nnkz*(nn-1)+1
             
tetra6(1,1)=nnkx*(nn-1)
tetra6(2,1)=nnky*(nn-1)
tetra6(3,1)=nnkz*(nn-1)
tetra6(1,2)=nnkx*(nn-1)
tetra6(2,2)=nnky*(nn-1)
tetra6(3,2)=nnkz*(nn-1)+1
tetra6(1,3)=nnkx*(nn-1)
tetra6(2,3)=nnky*(nn-1)+1
tetra6(3,3)=nnkz*(nn-1)+1
tetra6(1,4)=nnkx*(nn-1)+1
tetra6(2,4)=nnky*(nn-1)+1
tetra6(3,4)=nnkz*(nn-1)+1

if (cube(4,2-1)==1.and.cube(4,3-1)==1.and.cube(4,4-1)==1.and.&
   cube(4,5-1)==1.and.cube(4,6-1)==1.and.cube(4,7-1)==1.and.cube(4,8-1)==1.and.cube(4,9-1)==1) then
   tetra1=0
   tetra2=0
   tetra3=0
   tetra4=0
   tetra5=0
   tetra6=0
else
   if (cube(4,4-1)==1.and.cube(4,6-1)==1.and.cube(4,7-1)==1.and.cube(4,8-1)==1) then
      tetra2=0
      tetra3=0
      tetra4=0
      tetra5=0
      tetra6=0
   end if
   
   if ( cube(4,6-1)==1.and.cube(4,7-1)==1 ) then
      tetra2=0
      tetra5=0
      tetra6=0
   end if

   if ( cube(4,4-1)==1.and.cube(4,8-1)==1 ) then
      tetra3=0
      tetra4=0
      tetra6=0
   end if
   
   if (cube(4,5-1)==1.and.cube(4,7-1)==1.and.cube(4,8-1)==1.and.cube(4,9-1)==1) then
      tetra1=0
      tetra2=0
      tetra3=0
      tetra4=0
      tetra5=0
      tetra6=0
      tetra1(1,1)=nnkx*(nn-1)
      tetra1(2,1)=nnky*(nn-1)
      tetra1(3,1)=nnkz*(nn-1)
      tetra1(1,2)=nnkx*(nn-1)+1
      tetra1(2,2)=nnky*(nn-1)
      tetra1(3,2)=nnkz*(nn-1)
      tetra1(1,3)=nnkx*(nn-1)
      tetra1(2,3)=nnky*(nn-1)+1
      tetra1(3,3)=nnkz*(nn-1)
      tetra1(1,4)=nnkx*(nn-1)
      tetra1(2,4)=nnky*(nn-1)
      tetra1(3,4)=nnkz*(nn-1)+1
   end if
   
   if (cube(4,9-1)==1)  then
      tetra1=0
      tetra2=0
      tetra3=0
      tetra4=0
      tetra5=0
      tetra6=0
      tetra1(1,1)=nnkx*(nn-1)
      tetra1(2,1)=nnky*(nn-1)
      tetra1(3,1)=nnkz*(nn-1)
      tetra1(1,2)=nnkx*(nn-1)+1
      tetra1(2,2)=nnky*(nn-1)
      tetra1(3,2)=nnkz*(nn-1)
      tetra1(1,3)=nnkx*(nn-1)+1
      tetra1(2,3)=nnky*(nn-1)+1
      tetra1(3,3)=nnkz*(nn-1)
      tetra1(1,4)=nnkx*(nn-1)
      tetra1(2,4)=nnky*(nn-1)
      tetra1(3,4)=nnkz*(nn-1)+1
              
      tetra2(1,1)=nnkx*(nn-1)
      tetra2(2,1)=nnky*(nn-1)
      tetra2(3,1)=nnkz*(nn-1)
      tetra2(1,2)=nnkx*(nn-1)
      tetra2(2,2)=nnky*(nn-1)+1
      tetra2(3,2)=nnkz*(nn-1)
      tetra2(1,3)=nnkx*(nn-1)+1
      tetra2(2,3)=nnky*(nn-1)+1
      tetra2(3,3)=nnkz*(nn-1)
      tetra2(1,4)=nnkx*(nn-1)
      tetra2(2,4)=nnky*(nn-1)
      tetra2(3,4)=nnkz*(nn-1)+1
                 
      tetra3(1,1)=nnkx*(nn-1)+1
      tetra3(2,1)=nnky*(nn-1)
      tetra3(3,1)=nnkz*(nn-1)
      tetra3(1,2)=nnkx*(nn-1)+1
      tetra3(2,2)=nnky*(nn-1)+1
      tetra3(3,2)=nnkz*(nn-1)
      tetra3(1,3)=nnkx*(nn-1)
      tetra3(2,3)=nnky*(nn-1)
      tetra3(3,3)=nnkz*(nn-1)+1
      tetra3(1,4)=nnkx*(nn-1)+1
      tetra3(2,4)=nnky*(nn-1)
      tetra3(3,4)=nnkz*(nn-1)+1
                 
      tetra4(1,1)=nnkx*(nn-1)
      tetra4(2,1)=nnky*(nn-1)+1
      tetra4(3,1)=nnkz*(nn-1)
      tetra4(1,2)=nnkx*(nn-1)+1
      tetra4(2,2)=nnky*(nn-1)+1
      tetra4(3,2)=nnkz*(nn-1)
      tetra4(1,3)=nnkx*(nn-1)
      tetra4(2,3)=nnky*(nn-1)
      tetra4(3,3)=nnkz*(nn-1)+1
      tetra4(1,4)=nnkx*(nn-1)
      tetra4(2,4)=nnky*(nn-1)+1
      tetra4(3,4)=nnkz*(nn-1)+1
                 
      tetra5(1,1)=nnkx*(nn-1)+1
      tetra5(2,1)=nnky*(nn-1)+1
      tetra5(3,1)=nnkz*(nn-1)
      tetra5(1,2)=nnkx*(nn-1)
      tetra5(2,2)=nnky*(nn-1)
      tetra5(3,2)=nnkz*(nn-1)+1
      tetra5(1,3)=nnkx*(nn-1)+1
      tetra5(2,3)=nnky*(nn-1)
      tetra5(3,3)=nnkz*(nn-1)+1
      tetra5(1,4)=nnkx*(nn-1)
      tetra5(2,4)=nnky*(nn-1)+1
      tetra5(3,4)=nnkz*(nn-1)+1
!   write(10,*) tetra5(3,4)
   end if
end if
      
vol_tetra=0
ss=0
mat=0
roku=6.0d0


do ss=1,3
   mat(1,ss)=tetra1(1,ss+1)-tetra1(1,1)
   mat(2,ss)=tetra1(2,ss+1)-tetra1(2,1)
   mat(3,ss)=tetra1(3,ss+1)-tetra1(3,1)
end do
call matdet(mat,det)
vol_tetra(1)=det/roku

do ss=1,3
   mat(1,ss)=tetra2(1,ss+1)-tetra2(1,1)
   mat(2,ss)=tetra2(2,ss+1)-tetra2(2,1)
   mat(3,ss)=tetra2(3,ss+1)-tetra2(3,1)
end do
call matdet(mat,det)
vol_tetra(2)=det/roku

do ss=1,3
   mat(1,ss)=tetra3(1,ss+1)-tetra3(1,1)
   mat(2,ss)=tetra3(2,ss+1)-tetra3(2,1)
   mat(3,ss)=tetra3(3,ss+1)-tetra3(3,1)
end do
call matdet(mat,det)
vol_tetra(3)=det/roku

do ss=1,3
   mat(1,ss)=tetra4(1,ss+1)-tetra4(1,1)
   mat(2,ss)=tetra4(2,ss+1)-tetra4(2,1)
   mat(3,ss)=tetra4(3,ss+1)-tetra4(3,1)
end do
call matdet(mat,det)
vol_tetra(4)=det/roku

do ss=1,3
   mat(1,ss)=tetra5(1,ss+1)-tetra5(1,1)
   mat(2,ss)=tetra5(2,ss+1)-tetra5(2,1)
   mat(3,ss)=tetra5(3,ss+1)-tetra5(3,1)
end do
call matdet(mat,det)
vol_tetra(5)=det/roku

do ss=1,3
   mat(1,ss)=tetra6(1,ss+1)-tetra6(1,1)
   mat(2,ss)=tetra6(2,ss+1)-tetra6(2,1)
   mat(3,ss)=tetra6(3,ss+1)-tetra6(3,1)
end do
call matdet(mat,det)
vol_tetra(6)=det/roku

do i=1,6
   vol_tetra(i)=abs(vol_tetra(i))
end do

end subroutine tetra_apex_output



subroutine matdet(a,det)

implicit none

real(8),dimension(3,3),intent(in)::a
integer::i,j
real(8)::det

det=a(1,1)*a(2,2)*a(3,3)+a(3,2)*a(2,1)*a(1,3)+a(3,1)*a(1,2)*a(2,3)&
-a(3,1)*a(2,2)*a(1,3)-a(2,1)*a(1,2)*a(3,3)-a(3,2)*a(2,3)*a(1,1)
end subroutine matdet


subroutine tetra_energy_output(tetra1,tetra2,tetra3,tetra4,tetra5,tetra6,&
ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6)
real(8),dimension(3,4),intent(in)::tetra1,tetra2,tetra3,tetra4,tetra5,tetra6
integer::ss,indexx,i,j
real(8),dimension(100)::energy
real(8),dimension(4,100)::ene_tetra1,ene_tetra2,ene_tetra3,ene_tetra4,ene_tetra5,ene_tetra6
real(8),dimension(3,4)::a1,a2,a3,a4,a5,a6

ene_tetra1=0
ene_tetra2=0
ene_tetra3=0
ene_tetra4=0
ene_tetra5=0
ene_tetra6=0

do i=1,3
   do j=1,4
      a1(i,j)=tetra1(i,j)
      a2(i,j)=tetra2(i,j)
      a3(i,j)=tetra3(i,j)
      a4(i,j)=tetra4(i,j)
      a5(i,j)=tetra5(i,j)
      a6(i,j)=tetra6(i,j)
   end do
end do

energy=0
ene_tetra1=0
ene_tetra2=0
ene_tetra3=0
ene_tetra4=0
ene_tetra5=0
ene_tetra6=0

!$omp parallel 
!$omp do
do ss=1,4
   do indexx=1,6

call ham(a1(1,ss),a1(2,ss),a1(3,ss),energy)
      ene_tetra1(ss,indexx)=energy(indexx)*13.6d0
!      write(13,*) kx,ss,ene_tetra1(ss,indexx)

call ham(a2(1,ss),a2(2,ss),a2(3,ss),energy)
      ene_tetra2(ss,indexx)=energy(indexx)*13.6d0

call ham(a1(1,ss),a1(2,ss),a1(3,ss),energy)
      ene_tetra3(ss,indexx)=energy(indexx)*13.6d0

call ham(a1(1,ss),a1(2,ss),a1(3,ss),energy)
      ene_tetra4(ss,indexx)=energy(indexx)*13.6d0

call ham(a1(1,ss),a1(2,ss),a1(3,ss),energy)
      ene_tetra5(ss,indexx)=energy(indexx)*13.6d0

call ham(a1(1,ss),a1(2,ss),a1(3,ss),energy)
      ene_tetra6(ss,indexx)=energy(indexx)*13.6d0

   end do
end do
!$omp end do
!$omp end parallel

end subroutine tetra_energy_output

