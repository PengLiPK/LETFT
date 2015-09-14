! Calculate the true location from the inversion result
program migraloc

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: tmpt,tmpx,tmpy,tmpz,rx,ry,rz
real(kind=8) :: sx,sy,sz
real(kind=8) :: t0
integer :: i,ir,ns,nrcver,evnid
integer :: tpxnum,tpynum,updown
character(len=10) :: tmpstn,tmpname
character(len=70) :: inpf,outfile,tpfile




open(21,file='migraloc.inp',status='old')
read(21,*)inpf
read(21,*)outfile
read(21,*)ns
read(21,*)tpfile
read(21,*)tpminx,tpmaxx
read(21,*)tpminy,tpmaxy
read(21,*)tpxnum,tpynum
close(21)


! Read topography data
open(23,file=tpfile,status='old')
do i=1,tpxnum*tpynum
    read(23,*)topoxy(i)%x,topoxy(i)%y,topoz(i)
    topoxy(i)%num=i
end do
close(23)


! Calculate true location
open(30,file=inpf,status='old')
open(35,file=outfile,status='replace')
do i=1,ns
    read(30,*)nrcver,tmpname
    write(35,*)nrcver,tmpname
    read(30,*)tmpstn,sx,sy,sz
    write(35,2000)tmpstn,sx,sy,sz

    do ir=1,nrcver
        read(30,*)tmpt,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0

        ! Check if depth is above topography
        do while(.true.)
            call psurf(rx,ry,rz,topoxy,topoz,updown,&
                &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
            if(updown .eq. -1)then
                rz=rz+0.2
            else
                exit
            end if
        end do

        write(35,2001)tmpt,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0
    end do
    write(*,*)"Source",i,"is finished!"
end do
2000 format(a4,1x,f10.6,1x,f9.6,1x,f6.4)
2001 format(f10.7,1x,f10.6,1x,f9.6,1x,f6.4,1x,i4,1x,f10.6,1x,f9.6,1x,&
           &f7.4,1x,f10.7)


stop
end
