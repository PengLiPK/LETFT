!--------------------------------------------------------------------------
! This program convert the result of inversion to the 3D velocity model.
! Velocities of nodes with low DWS values keep the same.
!--------------------------------------------------------------------------
program cvtinvout

implicit none
integer,parameter :: imodel=50000
real(kind=8) :: nd(imodel,4)
real(kind=8) :: vel,tmpx,tmpy,tmpz
integer :: nx,ny,nz,velnum,outnodenum
integer :: totalnodenum,edgenodenum
integer :: im
character(len=70) :: oputmodel
character(len=70) :: invoutf
character(len=70) :: inpmodel

! Input parameters
!----------------------------------------------------

open(22,file='cvtinvout.inp',status='old')
read(22,*)invoutf
read(22,*)inpmodel
read(22,*)oputmodel
close(22)

!-----------------------------------------------------



! Read input model
!-----------------------------------------------------
open(32,file=inpmodel,status='old')
read(32,*)totalnodenum,nx,ny,nz
read(32,*)edgenodenum
do im=1,totalnodenum
    read(32,*)nd(im,1),nd(im,2),nd(im,3),nd(im,4)
end do
close(32)
!-----------------------------------------------------------


! Read invout file
!-----------------------------------------------------
open(35,file=invoutf,status='old')
read(35,*)outnodenum
read(35,*)
do im=1,outnodenum
    read(35,*)tmpx,tmpy,tmpz,vel,velnum
    nd(velnum,4)=vel
end do
close(35)
!-----------------------------------------------------


! Write output model
!-----------------------------------------------------
open(36,file=oputmodel,status='replace')
write(36,*)totalnodenum,nx,ny,nz
write(36,*)edgenodenum
do im=1,totalnodenum
    write(36,2000)nd(im,1),nd(im,2),nd(im,3),nd(im,4)
end do
2000 format(f10.6,1x,f10.6,1x,f9.6,1x,f8.6)
close(36)
!-----------------------------------------------------------


stop
end

