program cal_fdx

implicit none
real(kind=8) :: gval
real(kind=8) :: s(1000000)
real(kind=8) :: m(10000)
integer :: immax,iv
integer :: ifd,i,gid,gim,status1
character(len=70) :: FdFile
character(len=70) :: vfile
character(len=70) :: outf

open(21,file='cal_fdx.inp',status='old')
read(21,*)FdFile
read(21,*)vfile
read(21,*)outf
read(21,*)ifd,immax
close(21)
write(*,*)"ifd",ifd

open(23,file=vfile,status='old')
do iv=1,immax
    read(23,*)m(iv)
end do
close(23)



open(25,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=16)

s=0
i=0
do while(.true.)
    i=i+1
    read(25,rec=i,iostat=status1)gval,gid,gim
    if(status1/=0)exit
    write(26,*)gval,gid,gim
    s(gid)=s(gid)+gval*m(gim)
end do
close(25)

write(*,*)"ifd",ifd

open(27,file=outf,status='replace')
do i=1,ifd
    write(27,*)s(i)
end do
close(27)


stop
end
