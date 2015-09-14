! Calculate the true location from the inversion result
program dloc2loc_multi_notp

use strct
implicit none
real(kind=8) :: tmpt,tmpx,tmpy,tmpz,rx,ry,rz
real(kind=8) :: sx,sy,sz
real(kind=8) :: dloc(maxevn3d,3) 
real(kind=8) :: dt(maxevn3d),t0
real(kind=8) :: newloc(3)
real(kind=8) :: minx,maxx,miny,maxy,minz,maxz
real(kind=8) :: thrd,thrdt
integer :: i,ir,ns,nmodel,nrcver,evnid
character(len=10) :: tmpstn,tmpname
character(len=70) :: inpf,dlocfile,outfile




open(21,file='dloc2loc_multi_notp.inp',status='old')
read(21,*)inpf
read(21,*)dlocfile
read(21,*)outfile
read(21,*)ns,nmodel
read(21,*)thrd,thrdt
read(21,*)minx,maxx
read(21,*)miny,maxy
read(21,*)minz,maxz
close(21)


! Read inverted results of dlocation and dt0
open(26,file=dlocfile,status='old')
do i=1,nmodel
    read(26,*)dt(i)
    read(26,*)dloc(i,1)
    read(26,*)dloc(i,2)
    read(26,*)dloc(i,3)
    if(dt(i) .gt. thrdt)then
        dt(i)=thrdt
    else if(dt(i) .lt. -thrdt)then
        dt(i)=-thrdt
    end if
    do ir=1,3
        if(dloc(i,ir) .gt. thrd)then
            dloc(i,ir)=thrd
        else if(dloc(i,ir) .lt. -thrd)then
            dloc(i,ir)=-thrd
        end if
    end do
end do
close(26)


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

        call loc2ear(newloc(1),newloc(2),newloc(3),&
                    &rx,ry,rz,-dloc(evnid,1),&
                    &-dloc(evnid,2),-dloc(evnid,3))

        if(newloc(1) .lt. minx)then
            newloc(1)=minx
        else if(newloc(1) .gt. maxx)then
            newloc(1)=maxx
        end if
        if(newloc(2) .lt. miny)then
            newloc(2)=miny
        else if(newloc(2) .gt. maxy)then
            newloc(2)=maxy
        end if
        if(newloc(3) .lt. minz)then
            newloc(3)=minz
        else if(newloc(3) .gt. maxz)then
            newloc(3)=maxz
        end if

        write(35,2001)tmpt,tmpx,tmpy,tmpz,evnid,newloc(1),newloc(2),newloc(3),t0-dt(evnid)
    end do
end do
2000 format(a4,1x,f10.6,1x,f9.6,1x,f6.4)
2001 format(f10.7,1x,f10.6,1x,f9.6,1x,f6.4,1x,i4,1x,f10.6,1x,f9.6,1x,&
           &f7.4,1x,f10.7)


stop
end
