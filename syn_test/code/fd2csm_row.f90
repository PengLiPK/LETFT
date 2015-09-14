! Covert fd.dat to csm format with row number
program fd2csm_row


implicit none

real(kind=8) :: r8
integer :: irow(1000000)
integer :: ndata
integer :: irowtmp,icol,status1,i,j,k
character(len=70) :: fdfile
character(len=70) :: fdmtfile
character(len=70) :: csmfile

open(21,file='fd2csm_row.inp')
read(21,*)fdfile
read(21,*)fdmtfile
read(21,*)csmfile
read(21,*)ndata


open(103,file=fdfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(105,file=csmfile,status='replace',form='unformatted',&
    &access='direct',recl=12)

i=0
irow=0
k=1
do while(.true.)
    i=i+1
    read(103,rec=i,iostat=status1)r8,irowtmp,icol
    if(status1/=0)exit
    write(105,rec=i)r8,icol
    if(irowtmp .eq. k)then
        irow(k)=irow(k)+1
    else
        k=irowtmp
        irow(k)=irow(k)+1
    end if
end do


write(*,*)(i-1)*16

close(103)
close(105)

k=0
open(109,file=fdmtfile,status='replace')
do j=1,ndata
   k=k+irow(j)
   write(109,*)irow(j)
end do
write(*,*)k

stop
end
