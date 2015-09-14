! Initial input grid velocity model from 1D velocity.
! -----------------------------------------------------------------------
program initinpv

use strct
implicit none
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: v
integer :: nnode,edgenode
integer :: nlayer,vnl,nx,ny,nz
integer :: iv
character(len=70) :: outf
character(len=70) :: nodefile
character(len=70) :: vfile



open(21,file='initinpv.inp',status='old')
read(21,*)outf
read(21,*)nodefile
read(21,*)vfile


! Read 1D input velocity file
open(24,file=vfile,status='old')
read(24,*)nlayer
do iv=1,nlayer
    read(24,*)vinp(iv,1),vinp(iv,2)
end do


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


! Output grid velocity model.
open(30,file=outf,status='replace')
write(30,*)nnode,nx,ny,nz
write(30,*)edgenode
do iv=1,nnode
    call detlayer(vnl,nlayer,vinp(:,1),node(iv)%z)
    v=1.0d0/linterp(1.0d0/vinp(vnl,2),1.0d0/vinp(vnl+1,2),&
                   &vinp(vnl,1),vinp(vnl+1,1),node(iv)%z)
    write(30,2000)node(iv)%x,node(iv)%y,node(iv)%z,v
end do
2000 format(f10.5,1x,f10.5,1x,f10.5,1x,f10.6)
close(30)


stop
end



