program calmdres


implicit none
integer,parameter :: nmax=12500
real(kind=8) :: mfinal(nmax)
real(kind=8) :: dwsfinal(nmax),dws(nmax)
real(kind=8) :: gtmp
integer :: rctfinal(nmax),rct(nmax)
integer :: metagzero(nmax)
integer :: metag(nmax,2)
integer :: i,j,ir,status1
integer :: imt,im,idmax,immax,imtmax,izero,izeromax
integer :: rowstart,tmpmetag
character(len=70) :: gfile,metagf,outrctf,outdwsf


open(31,file='calmdres2.inp',status='old')
read(31,*)metagf
read(31,*)gfile
read(31,*)outrctf
read(31,*)outdwsf
read(31,*)idmax,immax
close(31)


! Read g
metag=0
open(40,file=metagf,status='old')
rowstart=1
imt=0
izero=0
do im=1,immax
    read(40,*)tmpmetag
    if(tmpmetag .ne. 0)then
        imt=imt+1
        metag(imt,1)=tmpmetag
        metag(imt,2)=rowstart
        rowstart=rowstart+metag(imt,1)
    else
        izero=izero+1
        metagzero(izero)=im
    end if
end do
imtmax=imt
izeromax=izero
close(40)


open(43,file=gfile,status='old',form='unformatted',&
    &access='direct',recl=12)
dws=0d0
rct=0
i=0
do imt=1,imtmax
    do j=1,metag(imt,1)
        i=i+1
        read(43,rec=i,iostat=status1)gtmp,ir
        if(status1/=0)exit
        dws(imt)=dws(imt)+abs(gtmp)
        if(gtmp .gt. 0d0)then
            rct(imt)=rct(imt)+1
        end if
    end do
end do


print *, 'Read g finished'
print *, 'immax=',immax,',  imtmax=',imtmax


print *, 'Calculate model finished'


! Calculate final ,rct,dws
mfinal=1d0
rctfinal=1
dwsfinal=1d0
do izero=1,izeromax
    mfinal(metagzero(izero))=0d0
    rctfinal(metagzero(izero))=0
    dwsfinal(metagzero(izero))=0d0
end do

imt=0
do im=1,immax
    if(mfinal(im) .ne. 0d0)then
        imt=imt+1
        rctfinal(im)=rct(imt)
        dwsfinal(im)=dws(imt)
    end if
end do


open(70,file=outrctf,status='replace')
do im=1,immax
    write(70,*)rctfinal(im)
end do
close(70)

open(75,file=outdwsf,status='replace')
do im=1,immax
    write(75,*)dwsfinal(im)
end do
close(75)


stop
end

