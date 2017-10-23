program test
implicit none
real(8)::c,det
c=4
call matdet(c,det)
write(*,*) det

end program test

subroutine matdet(c,d)
implicit none
real(8),intent(in)::c
real(8)::d
d=c
end subroutine matdet
