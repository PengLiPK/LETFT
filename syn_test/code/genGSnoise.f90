! This program generates random GS noises.
program genGSnoise

use strct
implicit none
integer,parameter :: maxns=10000
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: nsdata(maxns,3)
real(kind=8) :: mu,var,amp,ampt0
real(kind=8) :: tmpt,tmpx,tmpy,tmpz,rx,ry,rz,t0
real(kind=8) :: sx,sy,sz,noiseper
real(kind=8) :: minx,maxx,miny,maxy,minz,maxz
real(kind=8) :: rxtmp,rytmp,rztmp
integer :: evnid,is,ir,ns,nrcver
integer :: nsnum,model,para1
integer :: tpxnum,tpynum
integer :: updown
character(len=10) :: tmpstn,tmpname
character(len=70) :: inpf
character(len=70) :: outf
character(len=70) :: noisef
character(len=70) :: tpfile



open(21,file='genGSnoise.inp')
read(21,*)model
if(model .eq. 1)then
    read(21,*)outf
    read(21,*)mu,var,amp
    read(21,*)nsnum
else if(model .eq. 2)then
    read(21,*)inpf
    read(21,*)outf
    read(21,*)noisef
    read(21,*)ns
    read(21,*)minx,maxx
    read(21,*)miny,maxy
    read(21,*)minz,maxz
    read(21,*)mu,var,amp,ampt0
    read(21,*)para1
    read(21,*)tpfile
    read(21,*)tpminx,tpmaxx
    read(21,*)tpminy,tpmaxy
    read(21,*)tpxnum,tpynum

    ! Read topography data
    open(22,file=tpfile,status='old')
    do is=1,tpxnum*tpynum
        read(22,*)topoxy(is)%x,topoxy(is)%y,topoz(is)
        topoxy(is)%num=is
    end do
    close(22)
end if
close(21)


open(23,file=outf,status='replace')
if(model .eq. 1)then
    do is=1,nsnum
        call gsnoise(noiseper,mu,var)
        write(23,*)noiseper*amp
    end do
else if(model .eq. 2)then
    open(22,file=inpf,status='old')
    open(24,file=noisef,status='replace')

    ! Generate noise data
    do is=1,maxns
        do ir=1,3
            call gsnoise(noiseper,mu,var)
            nsdata(is,ir)=noiseper
        end do
    end do

    do is=1,ns
        read(22,*)nrcver,tmpname
        write(23,*)nrcver,tmpname
        read(22,*)tmpstn,sx,sy,sz
        write(23,2000)tmpstn,sx,sy,sz
        do ir=1,nrcver
            read(22,*)tmpt,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0
            if(para1 .eq. 1)then
                call gsnoise(noiseper,mu,var)
                tmpt=tmpt*(noiseper*amp+1.0d0)
            else if(para1 .eq. 2)then
                rxtmp=rx+nsdata(evnid,1)*amp*1.0d-2
                if(rxtmp .le. minx .or. rxtmp .ge. maxx)then
                    rx=rx-nsdata(evnid,1)*amp*1.0d-2
                else
                    rx=rxtmp
                end if
                rytmp=ry+nsdata(evnid,2)*amp*1.0d-2
                if(rytmp .le. miny .or. rytmp .ge. maxy)then
                    ry=ry-nsdata(evnid,2)*amp*1.0d-2
                else
                    ry=rytmp
                end if
                rztmp=rz+nsdata(evnid,3)*amp
                if(rztmp .le. minz .or. rztmp .ge. maxz)then
                    rz=rz-nsdata(evnid,3)*amp*1.0d-2
                else
                    rz=rztmp
                end if
                
                ! Avoid new depth is above topography.
                call psurf(rx,ry,rz,topoxy,topoz,updown,&
                    &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
                if(updown .eq. -1)then
                    rz=rz+abs(nsdata(evnid,3))*amp*1.0d-2
                end if

                t0=t0+nsdata(evnid,3)*ampt0
            end if
            write(23,2001)tmpt,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0
            write(24,*)noiseper*amp
        end do
    end do
2000 format(a4,1x,f9.5,1x,f8.5,1x,f5.3)
2001 format(f9.6,1x,f9.5,1x,f8.5,1x,f5.3,1x,i4,1x,f9.5,1x,f8.5,1x,&
                &f6.3,1x,f9.6)
    close(22)
    close(24)
end if
close(23)

stop
end
!---------------------------------------------------------------------




! Generate Gaussian noise by Central limit thorem method.
!---------------------------------------------------------------------
subroutine gsnoise(x,mu,var)

implicit none
real(kind=8) :: x,mu,var
real(kind=8),allocatable :: r(:)
integer,allocatable :: seed(:)
integer :: i,imax,clock,n


imax=1000
allocate(r(imax))


call random_seed(size=n)
allocate(seed(n))
call system_clock(count=clock)
seed=clock+(/(i-1,i=1,n)/)
call random_seed(put=seed)
deallocate(seed)
call random_number(r)

x=0d0
do i=1,imax
    x=x+r(i)
end do

! For uniform randoms in [0,1], mu=0.5, var=1/12
! Adjust x so mu=0, var=1
x=x-dble(imax)/2.0d0
x=x*sqrt(1.2d1/dble(imax))

! Adjust x to input mu and var
x=mu+sqrt(var)*x

return
end subroutine gsnoise
