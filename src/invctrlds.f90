
! Prevent large variation of the result in each iteration.
! -----------------------------------------------------------------------
program invctrlds

use strct
implicit none
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: sg(maxgrid3d)
real(kind=8) :: dws(maxgrid3d),dwsthrd
real(kind=8) :: vper
real(kind=8) :: tmps
real(kind=8) :: smin,smax
real(kind=8) :: vpermin,vpermax
integer :: nnode,edgenode
integer :: iv
integer :: nx,ny,nz
character(len=70) :: outf
character(len=70) :: inpf
character(len=70) :: nodefile
character(len=70) :: dwsfile



open(21,file='invctrlds.inp',status='old')
read(21,*)outf
read(21,*)inpf
read(21,*)nodefile
read(21,*)dwsfile
read(21,*)dwsthrd
read(21,*)smin,smax
read(21,*)vpermin,vpermax



! Read node file
open(25,file=nodefile,status='old')
read(25,*)nnode,nx,ny,nz
read(25,*)edgenode
do iv=1,nnode
    read(25,*)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t
    node(iv)%dxx=0d0
    node(iv)%num=iv
    node(iv)%stat=0
end do
close(25)


! Read s
open(26,file=inpf,status='old')
do iv=1,nnode
    read(26,*)sg(iv)
end do


! Read dws value
open(27,file=dwsfile,status='old')
do iv=1,nnode
    read(27,*)dws(iv)
end do


! Control slowness
open(30,file=outf,status='replace')
write(30,*)nnode,nx,ny,nz
write(30,*)edgenode
do iv=1,nnode
    if(dws(iv) .lt. dwsthrd)then
        write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t
    else
        tmps=1.0d0/node(iv)%t-sg(iv)
        if(tmps .gt. smax .or. tmps .lt. smin)then
            write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t
        else
            vper=((1.0d0/tmps)-node(iv)%t)/node(iv)%t
            if(vper .lt. vpermax .and. vper .gt. vpermin)then
                write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,1.0d0/tmps
            else if(vper .gt. vpermax)then
                write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t*(1.0d0+vpermax)
            else if(vper .lt. vpermin)then
                write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t*(1.0d0+vpermin)
            else
                write(*,*)"Error! vper,vpermax,vpermin:",vper,vpermax,vpermin
            end if
        end if
    end if
end do
2000 format(f10.5,1x,f10.5,1x,f10.5,1x,f10.6)
close(30)


stop
end



