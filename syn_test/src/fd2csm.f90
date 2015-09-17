! Covert fd.dat to csm format
program fd2csm


implicit none

real(kind=8) :: r8
integer :: icol(100000000)
integer :: ncol(1000000)
integer :: nmodel
integer :: irow,icoltmp,status1,i,j,k,m
character(len=70) :: fdfile
character(len=70) :: fdmtfile
character(len=70) :: csmfile

open(21,file='fd2csm.inp')
read(21,*)fdfile
read(21,*)fdmtfile
read(21,*)csmfile
read(21,*)nmodel


open(103,file=fdfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(105,file=csmfile,status='replace',form='unformatted',&
    &access='direct',recl=12)

i=0
do while(.true.)
    i=i+1
    read(103,rec=i,iostat=status1)r8,irow,icoltmp
    if(status1/=0)exit
    icol(i)=icoltmp
end do


write(*,*)(i-1)*16

m=0
ncol=0
do j=1,nmodel
    do k=1,i-1
        if(icol(k) .eq. j)then
            m=m+1
            read(103,rec=k,iostat=status1)r8,irow,icoltmp
            write(105,rec=m)r8,irow
            ncol(j)=ncol(j)+1
        end if
    end do
end do

close(103)
close(105)

k=0
open(109,file=fdmtfile,status='replace')
do j=1,nmodel
   k=k+ncol(j)
   write(109,*)ncol(j)
end do
write(*,*)k

stop
end
