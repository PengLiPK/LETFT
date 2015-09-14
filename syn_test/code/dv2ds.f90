
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program dv2ds

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: dvg(maxgrid3d)
real(kind=8) :: dsg(maxgrid3d)
real(kind=8) :: vg(maxgrid3d)
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: vtemp,vgtemp
real(kind=8) :: vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: topozbtm
real(kind=8) :: lon,lan,dpth
integer :: ngrd
integer :: nlayer
integer :: ix,iv,ig
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: outf
character(len=70) :: inpf
character(len=70) :: vfile
character(len=70) :: tpfile



open(21,file='dv2ds.inp',status='old')
read(21,*)outf
read(21,*)inpf
read(21,*)vfile
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



! Read ds,dv,v.
open(26,file=inpf,status='old')
open(30,file=outf,status='replace')
vg=0.0d0
read(26,*)ngrd
do ig=1,ngrd
    read(26,*)lon,lan,dpth,vgtemp
    call psurf(lon,lan,dpth,topoxy,topoz,updown,&
    &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
    if(updown .eq. -1)then
        dsg(ig)=0d0
        dvg(ig)=0d0
        call vfrom1d(vtemp,dpth,vinp,nlayer)
        vg(ig)=vtemp
    else if(updown .eq. 1)then
        dvg(ig)=vgtemp
        call vfrom1d(vtemp,dpth,vinp,nlayer)
        vg(ig)=vtemp+dvg(ig)
        dsg(ig)=(1.0d0/vg(ig))-(1.0d0/vtemp)
    else
        write(*,*)"error in Line113!"
    end if
    write(30,*)dsg(ig)
end do
close(26)
close(30)

stop
end
