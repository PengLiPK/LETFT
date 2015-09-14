! Split fd matrix into fdloc and fdvel parts.
program fdsplit



implicit none
integer,parameter :: nmax=100 
real(kind=8) :: f8(nmax)
integer :: i4(nmax)
integer :: nevn
integer :: i,j,k,status1
character(len=70) :: dfile,outlocdf,outveldf


open(21,file='fdsplit.inp',status='old')
read(21,*)dfile
read(21,*)outlocdf
read(21,*)outveldf
read(21,*)nevn

write(*,*)"Input file:",dfile
write(*,*)"Output file:",outlocdf,outveldf

open(103,file=dfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(105,file=outlocdf,status='replace',form='unformatted',&
    &access='direct',recl=16)
open(108,file=outveldf,status='replace',form='unformatted',&
    &access='direct',recl=16)
i=0
j=0
k=0
do while(.true.)
    i=i+1
    read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2)
    if(status1/=0)exit
    if(i4(2) .le. nevn*4)then
        j=j+1
        write(105,rec=j)f8(1),i4(1),i4(2)
    else
        k=k+1
        write(108,rec=k)f8(1),i4(1),i4(2)-nevn*4
    end if
end do
close(103)    
close(105)    
close(108)


stop
end
