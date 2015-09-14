
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program cutslice_reg2

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: dvg(maxgrid3d)
real(kind=8) :: dsg(maxgrid3d)
real(kind=8) :: vg(maxgrid3d)
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: vper,vtemp,vgtemp,sper
real(kind=8) :: dx,dy,dz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: x,y,z
real(kind=8) :: v,vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: tpx1,tpx2,tpy1,tpy2,tpsizex,tpsizey
real(kind=8) :: topozbtm,elev
integer :: nl(3),n(8)
integer :: nnode,edgenode
integer :: xnum,ynum,znum
integer :: nlayer
integer :: ix,iy,iz,iv,ig
integer :: itpx,itpy
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: outf
character(len=70) :: outpf
character(len=70) :: outvpf
character(len=70) :: inpf
character(len=70) :: nodefile
character(len=70) :: vfile
character(len=70) :: tpfile



open(21,file='cutslice_reg2b.inp',status='old')
read(21,*)outf
read(21,*)outpf
read(21,*)outvpf
read(21,*)inpf
read(21,*)nodefile
read(21,*)vfile
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

! Read 1D input velocity file
open(24,file=vfile,status='old')
read(24,*)nlayer
do iv=1,nlayer
    read(24,*)vinp(iv,1),vinp(iv,2),vinp(iv,3)
end do


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
            !call psurf(node(ig)%x,node(ig)%y,node(ig)%z,topoxy,topoz,updown,&
            !&tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
            !if(updown .eq. -1)then
            !    dsg(ig)=0d0
            !    dvg(ig)=0d0
            !    call vfrom1d(vtemp,node(ig)%z,vinp,nlayer)
            !    vg(ig)=vtemp
            !else if(updown .eq. 1)then
                dsg(ig)=vgtemp
                call vfrom1d(vtemp,node(ig)%z,vinp,nlayer)
                vg(ig)=vtemp/(-vtemp*dsg(ig)+1.0d0)
                dvg(ig)=(vg(ig)-vtemp)/vtemp
            !else
            !    write(*,*)"error in Line113!"
            !end if
        end do
    end do
end do
close(26)

xnum=nint((maxx-minx)/dx)+1
ynum=nint((maxy-miny)/dy)+1
znum=nint((maxz-minz)/dz)+1
open(30,file=outf,status='replace')
open(35,file=outpf,status='replace')
open(40,file=outvpf,status='replace')
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
               v=trilinear(vg(n(1)),vg(n(2)),vg(n(3)),vg(n(4)),&
                          &vg(n(5)),vg(n(6)),vg(n(7)),vg(n(8)),&
                          &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                          &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                          &x,y,z,1)
               !vper=trilinear(dvg(n(1)),dvg(n(2)),dvg(n(3)),dvg(n(4)),&
               !           &dvg(n(5)),dvg(n(6)),dvg(n(7)),dvg(n(8)),&
               !           &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
               !           &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
               !           &x,y,z,1)
               !sper=trilinear(dsg(n(1)),dsg(n(2)),dsg(n(3)),dsg(n(4)),&
               !           &dsg(n(5)),dsg(n(6)),dsg(n(7)),dsg(n(8)),&
               !           &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
               !           &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
               !           &x,y,z,1)
               vper=dvg(n(8))
               sper=dsg(n(8))
               write(30,*)x,y,z,v
               write(35,*)x,y,z,sper
               write(40,*)x,y,z,vper
           else
               call psurf(x,y,z,topoxy,topoz,updown,&
               &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
               if(updown .eq. -1)then
                   write(30,*)x,y,z," NaN"
                   write(35,*)x,y,z," NaN"
                   write(40,*)x,y,z," NaN"
               else if(updown .eq. 1)then
                   v=trilinear(vg(n(1)),vg(n(2)),vg(n(3)),vg(n(4)),&
                              &vg(n(5)),vg(n(6)),vg(n(7)),vg(n(8)),&
                              &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                              &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                              &x,y,z,1)
                   !vper=trilinear(dvg(n(1)),dvg(n(2)),dvg(n(3)),dvg(n(4)),&
                   !           &dvg(n(5)),dvg(n(6)),dvg(n(7)),dvg(n(8)),&
                   !           &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                   !           &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                   !           &x,y,z,1)
                   !sper=trilinear(dsg(n(1)),dsg(n(2)),dsg(n(3)),dsg(n(4)),&
                   !           &dsg(n(5)),dsg(n(6)),dsg(n(7)),dsg(n(8)),&
                   !           &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                   !           &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                   !           &x,y,z,1)
                   vper=dvg(n(8))
                   sper=dsg(n(8))
                   write(30,*)x,y,z,v
                   write(35,*)x,y,z,sper
                   write(40,*)x,y,z,vper
               else
                   write(*,*)"Paramter updown(refined nodes) is not -1,1.&
                             & Error in subroutine psurf in strct.f90."
               end if
           end if

        end do
    end do
end do
close(30)
close(35)
close(40)

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



