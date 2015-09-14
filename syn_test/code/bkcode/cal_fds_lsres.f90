program cal_fds_lsres

implicit none
real(kind=8) :: gval,sums,sumsper,sper,sumlsres
real(kind=8) :: s(1000000)
real(kind=8) :: t(1000000),t0(1000000)
real(kind=8) :: lsres(10000)
real(kind=8) :: node(10000)
integer :: nnode,edgenode,iv
integer :: ifd,i,gid,gim,status1
character(len=70) :: FdFile
character(len=70) :: vfile
character(len=70) :: dfile,d0file
character(len=70) :: outf
character(len=70) :: outflsres

open(21,file='cal_fds_lsres.inp',status='old')
read(21,*)FdFile
read(21,*)vfile
read(21,*)dfile
read(21,*)d0file
read(21,*)outf
read(21,*)outflsres
read(21,*)ifd,nnode
close(21)
write(*,*)"ifd",ifd


open(22,file=dfile,status='old')
do i=1,ifd
    read(22,*)t(i)
end do
close(22)

open(23,file=d0file,status='old')
do i=1,ifd
    read(23,*)t0(i)
end do
close(23)


open(24,file=vfile,status='old')
do iv=1,nnode
    read(24,*)node(iv)
end do
close(24)



open(25,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=16)

open(26,file='fdasc.txt',status='replace')
s=0
i=0

! Calculate G0*dm
do while(.true.)
    i=i+1
    read(25,rec=i,iostat=status1)gval,gid,gim
    if(status1/=0)exit
    write(26,*)gval,gid,gim
    !if(gim .gt. 48)then
        s(gid)=s(gid)+gval*node(gim)
    !end if
end do

!lsres=0
!i=0
!do while(.true.)
!    i=i+1
!    read(25,rec=i,iostat=status1)gval,gid,gim
!    if(status1/=0)exit
!    if(gim .gt. 48)then
!        lsres(gim)=lsres(gim)+gval*s(gid)
!    end if
!end do

close(25)
close(26)

write(*,*)"ifd",ifd

open(27,file=outf,status='replace')
sums=0
sumsper=0
do i=1,ifd
    s(i)=s(i)+t0(i)-t(i)
    sums=sums+abs(s(i))
    sper=s(i)/t(i)
    sumsper=sumsper+abs(sper)
    write(27,*)sper,s(i)
end do
sums=sums/ifd
sumsper=sumsper/ifd
write(27,*)"ave_resper=:",sumsper
write(27,*)"ave_res=:",sums
close(27)

!open(28,file=outflsres,status='replace')
!sumlsres=0
!do i=1,nnode
!    write(28,*)lsres(i)
!    sumlsres=sumlsres+abs(lsres(i))
!end do
!sumlsres=sumlsres/nnode
!write(28,*)"ave_lsres=:",sumlsres
!close(28)

stop
end
