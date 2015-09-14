! This module defines some struct type variables which will be used in 
! main program and subroutine.
!------------------------------------------------------------------------

module mdcutprb4inv
    implicit none
    integer,parameter :: imodel=50000

    type gstrct
        real(kind=8) :: gen
        integer :: im
        integer :: id
    end type
end module mdcutprb4inv
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
! This program cut nodes with small dws values.
! The outputs nodefile and fd file will be used in inversion.
!--------------------------------------------------------------------------
program cutpro4inv


use mdcutprb4inv
implicit none
type(gstrct) :: g
real(kind=8) :: dws(imodel)
real(kind=8) :: nd(imodel,4),newnd(imodel,4)
real(kind=8) :: dwsthrd
integer :: ndnum(imodel),imnum(imodel)
integer :: totalnodenum,edgenodenum
integer :: inm,im
integer :: ifd,ifdnew
integer :: status1
character(len=70) :: dwsfile
character(len=70) :: fdfile
character(len=70) :: oputmodel
character(len=70) :: oputfd
character(len=70) :: inpmodel

! Input parameters
!----------------------------------------------------

open(22,file='cutpro4inv.inp')
read(22,*)fdfile
read(22,*)dwsfile
read(22,*)inpmodel
read(22,*)oputfd
read(22,*)oputmodel
read(22,*)dwsthrd
close(22)

!-----------------------------------------------------


! Read input model
!-----------------------------------------------------
open(32,file=inpmodel,status='old')
read(32,*)totalnodenum
read(32,*)edgenodenum
do im=1,totalnodenum
    read(32,*)nd(im,1),nd(im,2),nd(im,3),nd(im,4)
end do
close(32)
!-----------------------------------------------------------


! Read dws values
!-----------------------------------------------------
open(35,file=dwsfile,status='old')
do im=1,totalnodenum
    read(35,*)dws(im)
end do
close(35)
!-----------------------------------------------------------


! Write node with dws > dwsthrd.
!-----------------------------------------------------
ndnum=0
imnum=0
inm=0
do im=1,totalnodenum
    if(dws(im) .gt. dwsthrd)then
        inm=inm+1
        newnd(inm,:)=nd(im,:)        
        ndnum(inm)=im
        imnum(im)=inm
    end if
end do

open(37,file=oputmodel,status='replace')
write(37,*)inm
im=0
write(37,*)im
do im=1,inm
    write(37,2000)newnd(im,1),newnd(im,2),newnd(im,3),newnd(im,4),ndnum(im)
end do
2000 format(f10.6,1x,f10.6,1x,f9.6,1x,f8.6,1x,i6)
close(37)
!-----------------------------------------------------------


! Reading G and write new G
!-----------------------------------------------------------
open(41,file=fdfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(42,file=oputfd,status='replace',form='unformatted',&
    &access='direct',recl=16)
ifd=0
ifdnew=0
do while(.true.)
    ifd=ifd+1
    read(41,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    if(dws(g%im) .gt. dwsthrd)then
        ifdnew=ifdnew+1
        write(42,rec=ifdnew)g%gen,g%id,imnum(g%im)
    end if
end do
close(41)
close(42)

end

