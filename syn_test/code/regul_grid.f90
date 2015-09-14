
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program regul_grid

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: vtemp,vtemp1
real(kind=8) :: minx(maxgrd1d),maxx(maxgrd1d)
real(kind=8) :: miny(maxgrd1d),maxy(maxgrd1d)
real(kind=8) :: minz(maxgrd1d),maxz(maxgrd1d)
real(kind=8) :: gsizex(maxgrd1d),gsizey(maxgrd1d),gsizez(maxgrd1d)
real(kind=8) :: x(maxgrd1d),y(maxgrd1d),z(maxgrd1d)
real(kind=8) :: vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: topozbtm
integer :: nxb,nyb,nzb
integer :: xbgnum,ybgnum,zbgnum
integer :: xgnum,ygnum,zgnum
integer :: nlayer
integer :: ix,iy,iz,iv
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: inpf
character(len=70) :: outf
character(len=70) :: outf1
character(len=70) :: vfile
character(len=70) :: tpfile



open(21,file='regul_grid.inp',status='old')
read(21,*)inpf
read(21,*)outf
read(21,*)outf1
read(21,*)vfile
read(21,*)tpfile
read(21,*)tpminx,tpmaxx
read(21,*)tpminy,tpmaxy
read(21,*)tpxnum,tpynum
read(21,*)topozbtm
read(21,*)vair


! Read the distribute of nodes.
open(22,file=inpf,status='old')
read(22,*)nxb
do ix=1,nxb
    read(22,*)minx(ix),maxx(ix),gsizex(ix)
end do

read(22,*)nyb
do iy=1,nyb
    read(22,*)miny(iy),maxy(iy),gsizey(iy)
end do

read(22,*)nzb
do iz=1,nzb
    read(22,*)minz(iz),maxz(iz),gsizez(iz)
end do
close(22)

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


! Initial x,y,z
xgnum=0
do ix=1,nxb
    xbgnum=nint((maxx(ix)-minx(ix))/gsizex(ix))+1
    do iv=1,xbgnum
        xgnum=xgnum+1
        x(xgnum)=minx(ix)+dble(iv-1)*gsizex(ix)
    end do
end do

ygnum=0
do iy=1,nyb
    ybgnum=nint((maxy(iy)-miny(iy))/gsizey(iy))+1
    do iv=1,ybgnum
        ygnum=ygnum+1
        y(ygnum)=miny(iy)+dble(iv-1)*gsizey(iy)
    end do
end do

zgnum=0
do iz=1,nzb
    zbgnum=nint((maxz(iz)-minz(iz))/gsizez(iz))+1
    do iv=1,zbgnum
        zgnum=zgnum+1
        z(zgnum)=minz(iz)+dble(iv-1)*gsizez(iz)
    end do
end do


! Generate grid file
open(30,file=outf,status='replace')
open(31,file=outf1,status='replace')
write(30,*)xgnum*ygnum*zgnum,xgnum,ygnum,zgnum
write(30,*)0
write(31,*)xgnum*ygnum*zgnum,xgnum,ygnum,zgnum
write(31,*)0
do iz=1,zgnum
    do iy=1,ygnum
        do ix=1,xgnum
            if(z(iz) .gt. topozbtm)then
                call vfrom1d(vtemp,z(iz),vinp,nlayer)
            else
                call psurf(x(ix),y(iy),z(iz),topoxy,topoz,updown,&
                &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
                if(updown .eq. -1)then
                    vtemp=vair
                else if(updown .eq. 1)then
                    call vfrom1d(vtemp,z(iz),vinp,nlayer)
                else
                    write(*,*)"Paramter updown(refined nodes) is not -1,1.&
                              & Error in subroutine psurf in strct.f90."
                end if
            end if
            call vfrom1d(vtemp1,z(iz),vinp,nlayer)
            write(31,*)x(ix),y(iy),z(iz),vtemp1
            write(30,*)x(ix),y(iy),z(iz),vtemp
        end do
    end do
end do
close(30)
close(31)

stop
end
