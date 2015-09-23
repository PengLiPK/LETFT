
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program cut_ressl_reg2

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: dsg(maxgrid3d)
real(kind=8) :: vg(maxgrid3d)
real(kind=8) :: vgtemp,sper
real(kind=8) :: dx,dy,dz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: x,y,z
real(kind=8) :: vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: topozbtm
integer :: nl(3),n(8)
integer :: nnode,edgenode
integer :: xnum,ynum,znum
integer :: ix,iy,iz,iv,ig
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: outpf
character(len=70) :: inpf
character(len=70) :: nodefile
character(len=70) :: tpfile



open(21,file='cut_ressl_reg2.inp',status='old')
read(21,*)outpf
read(21,*)inpf
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
open(25,file=nodefile,status='old')
read(25,*)nnode,nl(1),nl(2),nl(3)
read(25,*)edgenode
do iv=1,nnode
    read(25,*)node(iv)%x,node(iv)%y,node(iv)%z,node(iv)%t
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




! Read ds,dv,v.
open(26,file=inpf,status='old')
vg=0.0d0
write(*,*)nl
do iz=1,nl(3)
    do iy=1,nl(2)
        do ix=1,nl(1)
            ig=nl(2)*nl(1)*(iz-1)+nl(1)*(iy-1)+ix
            read(26,*)vgtemp
            dsg(ig)=vgtemp
        end do
    end do
end do
close(26)

xnum=nint((maxx-minx)/dx)+1
ynum=nint((maxy-miny)/dy)+1
znum=nint((maxz-minz)/dz)+1
open(35,file=outpf,status='replace')
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
               sper=trilinear2(dsg(n(1)),dsg(n(2)),dsg(n(3)),dsg(n(4)),&
                          &dsg(n(5)),dsg(n(6)),dsg(n(7)),dsg(n(8)),&
                          &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                          &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                          &x,y,z,topoxy,topoz,tpminx,tpmaxx,tpminy,&
                          &tpmaxy,tpxnum,tpynum,1d5)
               write(35,*)x,y,z,sper
           else
               call psurf(x,y,z,topoxy,topoz,updown,&
               &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
               if(updown .eq. -1)then
                   write(35,*)x,y,z," NaN"
               else if(updown .eq. 1)then
                   sper=trilinear2(dsg(n(1)),dsg(n(2)),dsg(n(3)),dsg(n(4)),&
                              &dsg(n(5)),dsg(n(6)),dsg(n(7)),dsg(n(8)),&
                              &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                              &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                              &x,y,z,topoxy,topoz,tpminx,tpmaxx,tpminy,&
                              &tpmaxy,tpxnum,tpynum,1d5)
                   write(35,*)x,y,z,sper
               else
                   write(*,*)"Paramter updown(refined nodes) is not -1,1.&
                             & Error in subroutine psurf in strct.f90."
               end if
           end if

        end do
    end do
end do
close(35)


stop
end



