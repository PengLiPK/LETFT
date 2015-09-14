
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program cutvdiff_reg1

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: sg(maxgrid3d)
real(kind=8) :: dx,dy,dz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: x,y,z
real(kind=8) :: s,v,vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: tpx1,tpx2,tpy1,tpy2,tpsizex,tpsizey
real(kind=8) :: topozbtm,elev
integer :: nl(3),n(8)
integer :: nnode,edgenode
integer :: xnum,ynum,znum
integer :: ix,iy,iz,iv
integer :: itpx,itpy
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: outf
character(len=70) :: nodefile
character(len=70) :: tpfile



open(21,file='cutvdiff_reg1.inp',status='old')
read(21,*)outf
read(21,*)nodefile
read(21,*)minx,maxx
read(21,*)miny,maxy
read(21,*)minz,maxz
read(21,*)dx,dy,dz
read(21,*)tpfile
read(21,*)tpminx,tpmaxx
read(21,*)tpminy,tpmaxy
read(21,*)tpxnum,tpynum
read(21,*)topozbtm
read(21,*)vair

! Read topography data
open(23,file=tpfile,status='old')
do ix=1,tpxnum*tpynum
    read(23,*)topoxy(ix)%x,topoxy(ix)%y,topoz(ix)
    topoxy(ix)%num=ix
end do
close(23)


! Read node file
sg=0d0
open(25,file=nodefile,status='old')
read(25,*)nnode,nl(1),nl(2),nl(3)
read(25,*)edgenode
do iv=1,nnode
    read(25,*)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t
    call psurf(node(iv)%x,node(iv)%y,node(iv)%z,topoxy,topoz,updown,&
        &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
    if(updown .eq. -1)then
        sg(iv)=0d0
    else
        sg(iv)=node(iv)%t
    end if
    node(iv)%dxx=0d0
    node(iv)%num=iv
    node(iv)%stat=0
end do
close(25)

do iv=1,nl(1)
    ix=iv
    layer(iv,1)=node(ix)%x
end do

do iv=1,nl(2)
    ix=(iv-1)*nl(1)+1
    layer(iv,2)=node(ix)%y
end do

do iv=1,nl(3)
    ix=(iv-1)*nl(1)*nl(2)+1
    layer(iv,3)=node(ix)%z
end do



xnum=nint((maxx-minx)/dx)+1
ynum=nint((maxy-miny)/dy)+1
znum=nint((maxz-minz)/dz)+1
open(30,file=outf,status='replace')
write(*,*)xnum,ynum,znum

do iz=1,znum
    do iy=1,ynum
        do ix=1,xnum
           x=minx+dble(ix-1)*dx
           y=miny+dble(iy-1)*dy 
           z=minz+dble(iz-1)*dz
           call locatcood3d2(n,nl(1),layer(1:maxgrd1d,1),nl(2),layer(1:maxgrd1d,2)&
                            &,nl(3),layer(1:maxgrd1d,3),x,y,z)
           if(z .gt. topozbtm)then
               s=trilinear(sg(n(1)),sg(n(2)),sg(n(3)),sg(n(4)),&
                          &sg(n(5)),sg(n(6)),sg(n(7)),sg(n(8)),&
                          &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                          &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                          &x,y,z,1)
               v=s
               write(30,*)x,y,z,v
           else
               call psurf(x,y,z,topoxy,topoz,updown,&
               &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
               if(updown .eq. -1)then
                   write(30,*)x,y,z," NaN"
               else if(updown .eq. 1)then
                   s=trilinear(sg(n(1)),sg(n(2)),sg(n(3)),sg(n(4)),&
                              &sg(n(5)),sg(n(6)),sg(n(7)),sg(n(8)),&
                              &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                              &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                              &x,y,z,1)
                   v=s
                   write(30,*)x,y,z,v
               else
                   write(*,*)"Paramter updown(refined nodes) is not -1,1.&
                             & Error in subroutine psurf in strct.f90."
               end if
           end if

        end do
    end do
end do
close(30)

! Calculate elevation lines
if(znum .ne. 1)then
    tpsizex=(tpmaxx-tpminx)/dble(tpxnum-1)
    tpsizey=(tpmaxy-tpminy)/dble(tpynum-1)
    open(45,file='elev_line.txt',status='replace')
    do iy=1,ynum
        do ix=1,xnum
           x=minx+dble(ix-1)*dx
           y=miny+dble(iy-1)*dy 
           itpx=int((x-tpminx)/tpsizex)+1
           itpy=int((y-tpminy)/tpsizey)+1
           tpx1=tpminx+dble(itpx-1)*tpsizex
           tpy1=tpminy+dble(itpy-1)*tpsizey
           tpx2=tpminx+dble(itpx)*tpsizex
           tpy2=tpminy+dble(itpx)*tpsizey
           elev=bilinear(topoz(itpx+(itpy-1)*tpxnum),&
                        &topoz(itpx+1+(itpy-1)*tpxnum),&
                        &topoz(itpx+itpy*tpxnum),&
                        &topoz(itpx+1+itpy*tpxnum),&
                        &tpx1,tpy1,tpx2,tpy2,x,y)
           write(45,*)x,y,elev
        end do
    end do
    close(45)
end if


stop
end



