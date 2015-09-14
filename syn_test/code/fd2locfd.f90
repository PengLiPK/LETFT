! Convert fd.dat to fdloc.dat.
! npr = 1: fdloc.dat for inverting x,y,z,t0
!       2: fdloc.dat for inverting x,y,z
program fd2locfd

implicit none
integer,parameter :: nmax=100 
real(kind=8) :: f8(nmax)
integer :: i4(nmax),npr,nevn
integer :: i,j,status1
character(len=70) :: dfile,outf


open(21,file='fd2locfd.inp',status='old')
read(21,*)dfile
read(21,*)outf
read(21,*)nevn
read(21,*)npr

write(*,*)"Input file:",dfile
write(*,*)"Output file:",outf

open(103,file=dfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(108,file=outf,status='replace',form='unformatted',&
    &access='direct',recl=16)
i=0
j=0
if(npr .eq. 1)then
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2)
        if(status1/=0)exit
        if(i4(2) .le. 4*nevn .and. i4(2) .ge. 1)then
            j=j+1
            write(108,rec=j)f8(1),i4(1),i4(2)
        end if
    end do
else if(npr .eq. 2)then
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2)
        if(status1/=0)exit
        if(i4(2) .le. 4*nevn .and. mod(i4(2),4) .ne. 1)then
            j=j+1
            write(108,rec=j)f8(1),i4(1),i4(2)-1
        end if
    end do
end if
    
close(103)    
close(108)


stop
end
