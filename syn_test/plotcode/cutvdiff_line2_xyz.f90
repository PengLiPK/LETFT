
! This program caculate velocity in each grid.
! -----------------------------------------------------------------------
program cutvdiff_line2_xyz

use strct
implicit none
type(tstrct) :: topoxy(maxgrid)
type(tstrct3d) :: node(maxgrid3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: sg(maxgrid3d)
real(kind=8) :: mdres(maxgrid3d)
real(kind=8) :: lon(maxvel),lan(maxvel),llen(maxvel)
real(kind=8) :: xlmax,ylmax,zlmax,lontmp,lantmp,ztmp
real(kind=8) :: ndtmpx,ndtmpy,ndtmpz
real(kind=8) :: ndprj,nddist,ndthrsh
real(kind=8) :: dl,dz,tmplen
real(kind=8) :: minz,maxz
real(kind=8) :: z,xx,yy,hrlen
real(kind=8) :: s,v,vair
real(kind=8) :: tpminx,tpmaxx
real(kind=8) :: tpminy,tpmaxy
real(kind=8) :: tpx1,tpx2,tpy1,tpy2,tpsizex,tpsizey
real(kind=8) :: topozbtm,elev,thrshd,resl
integer :: nl(3),n(8)
integer :: nnode,edgenode
integer :: npts,lnum,znum
integer :: il,ihr,iz,iv
integer :: itpx,itpy
integer :: tpxnum,tpynum
integer :: updown
character(len=70) :: outf,outresf,outndprjf
character(len=70) :: nodefile,mdresf,lfile
character(len=70) :: tpfile



open(21,file='cutvdiff_line2_xyz.inp',status='old')
read(21,*)outf
read(21,*)outresf
read(21,*)outndprjf,ndthrsh
read(21,*)nodefile
read(21,*)mdresf,thrshd
read(21,*)lfile
read(21,*)minz,maxz
read(21,*)dl,dz
read(21,*)tpfile
read(21,*)tpminx,tpmaxx
read(21,*)tpminy,tpmaxy
read(21,*)tpxnum,tpynum
read(21,*)topozbtm
read(21,*)vair

! Read topography data
open(23,file=tpfile,status='old')
do iv=1,tpxnum*tpynum
    read(23,*)topoxy(iv)%x,topoxy(iv)%y,topoz(iv)
    topoxy(iv)%num=iv
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
        sg(iv)=vair
        !sg(iv)=0d0
    else
        sg(iv)=node(iv)%t
    end if
    node(iv)%dxx=0d0
    node(iv)%num=iv
    node(iv)%stat=0
end do
close(25)

do iv=1,nl(1)
    ihr=iv
    layer(iv,1)=node(ihr)%x
end do

do iv=1,nl(2)
    ihr=(iv-1)*nl(1)+1
    layer(iv,2)=node(ihr)%y
end do

do iv=1,nl(3)
    ihr=(iv-1)*nl(1)*nl(2)+1
    layer(iv,3)=node(ihr)%z
end do


! Read model res file
open(27,file=mdresf,status='old')
do iv=1,nnode
    read(27,*)mdres(iv)
end do
close(27)

! Read line file
open(28,file=lfile,status='old')
read(28,*)npts
do iv=1,npts
    read(28,*)lon(iv),lan(iv)
end do
close(28)

znum=nint((maxz-minz)/dz)+1
open(30,file=outf,status='replace')
open(35,file=outresf,status='replace')
open(38,file=outndprjf,status='replace')

tmplen=0d0
do il=1,npts-1
    call ear2loc(xlmax,ylmax,zlmax,lon(il+1),lan(il+1),0d0,lon(il),lan(il),0d0)
    llen(il)=sqrt(xlmax*xlmax+ylmax*ylmax)
    lnum=int(llen(il)/dl)
    if(il .gt. 1)then
        tmplen=tmplen+llen(il-1)
    end if
    do ihr=1,lnum
        xx=dl*dble(ihr-1)*xlmax/llen(il)
        yy=dl*dble(ihr-1)*ylmax/llen(il)
        call loc2ear(lontmp,lantmp,ztmp,lon(il),lan(il),0d0,xx,yy,0d0)
        hrlen=dl*dble(ihr-1)+tmplen
        do iz=1,znum
           z=minz+dble(iz-1)*dz
           call locatcood3d2(n,nl(1),layer(1:maxgrd1d,1),nl(2),layer(1:maxgrd1d,2)&
                            &,nl(3),layer(1:maxgrd1d,3),lontmp,lantmp,z)
           resl=trilinear2(mdres(n(1)),mdres(n(2)),mdres(n(3)),mdres(n(4)),&
                      &mdres(n(5)),mdres(n(6)),mdres(n(7)),mdres(n(8)),&
                      &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                      &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                      &lontmp,lantmp,z,topoxy,topoz,tpminx,tpmaxx,tpminy,&
                      &tpmaxy,tpxnum,tpynum,1d5)
           if(resl .ge. thrshd)then
               if(z .gt. topozbtm)then
                   s=trilinear(sg(n(1)),sg(n(2)),sg(n(3)),sg(n(4)),&
                              &sg(n(5)),sg(n(6)),sg(n(7)),sg(n(8)),&
                              &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                              &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                              &lontmp,lantmp,z,1)
                   v=s
                   write(30,*)hrlen,z,v
                   write(35,*)hrlen,z,resl
               else
                   call psurf(lontmp,lantmp,z,topoxy,topoz,updown,&
                   &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
                   if(updown .eq. -1)then
                       write(30,*)hrlen,z," NaN"
                       write(35,*)hrlen,z," NaN"
                   else if(updown .eq. 1)then
                       s=trilinear(sg(n(1)),sg(n(2)),sg(n(3)),sg(n(4)),&
                                  &sg(n(5)),sg(n(6)),sg(n(7)),sg(n(8)),&
                                  &node(n(1))%x,node(n(1))%y,node(n(1))%z,&
                                  &node(n(8))%x,node(n(8))%y,node(n(8))%z,&
                                  &lontmp,lantmp,z,1)
                       v=s
                       write(30,*)hrlen,z,v
                       write(35,*)hrlen,z,resl
                   else
                       write(*,*)"Paramter updown(refined nodes) is not -1,1.&
                                 & Error in subroutine psurf in strct.f90."
                   end if
               end if
           else
               write(30,*)hrlen,z," NaN"
               write(35,*)hrlen,z,resl
           end if

        end do
    end do

    ! Projecting nodes closed to profile
    do iv=1,nnode
        !if(node(iv)%x .ge. min(lon(il),lon(il+1)) .and. &
        !  &node(iv)%x .le. max(lon(il),lon(il+1)) .and. &
        !  &node(iv)%y .ge. min(lan(il),lan(il+1)) .and. &
        !  &node(iv)%y .le. max(lan(il),lan(il+1)))then
            
            call psurf(node(iv)%x,node(iv)%y,node(iv)%z,topoxy,topoz,updown,&
            &tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)
            if(updown .eq. 1)then
                call ear2loc(ndtmpx,ndtmpy,ndtmpz,node(iv)%x,node(iv)%y,&
                            &0d0,lon(il),lan(il),0d0)

                ! Dot product A*B=|A|*|B|*cos(c)
                ndprj=(ndtmpx*xlmax+ndtmpy*ylmax)/llen(il)
                nddist=sqrt(ndtmpx*ndtmpx+ndtmpy*ndtmpy-ndprj*ndprj)

                if(nddist .le. ndthrsh .and. ndprj .ge. 0d0 .and. &
                  &ndprj .le. llen(il))then
                    write(38,*)ndprj+tmplen,node(iv)%z
                end if
            end if
        !end if
    end do
end do

tmplen=0d0
do il=1,npts-1
    tmplen=tmplen+llen(il)
    write(30,*)tmplen
end do

close(30)
close(35)
close(38)

! Calculate elevation lines
tpsizex=(tpmaxx-tpminx)/dble(tpxnum-1)
tpsizey=(tpmaxy-tpminy)/dble(tpynum-1)
open(45,file='elev_line.txt',status='replace')
open(50,file='elev_poly.txt',status='replace')
write(50,*)0.000,0.000
tmplen=0d0
do il=1,npts-1
    call ear2loc(xlmax,ylmax,zlmax,lon(il+1),lan(il+1),0d0,lon(il),lan(il),0d0)
    llen(il)=sqrt(xlmax**2+ylmax**2)
    lnum=int(llen(il)/dl)
    if(il .gt. 1)then
        tmplen=tmplen+llen(il-1)
    end if
    do ihr=1,lnum
        xx=dl*dble(ihr-1)*xlmax/llen(il)
        yy=dl*dble(ihr-1)*ylmax/llen(il)
        call loc2ear(lontmp,lantmp,ztmp,lon(il),lan(il),0d0,xx,yy,0d0)
        hrlen=dl*dble(ihr-1)+tmplen
        itpx=int((lontmp-tpminx)/tpsizex)+1
        itpy=int((lantmp-tpminy)/tpsizey)+1
        tpx1=tpminx+dble(itpx-1)*tpsizex
        tpy1=tpminy+dble(itpy-1)*tpsizey
        tpx2=tpminx+dble(itpx)*tpsizex
        tpy2=tpminy+dble(itpx)*tpsizey
        elev=bilinear(topoz(itpx+(itpy-1)*tpxnum),&
                    &topoz(itpx+1+(itpy-1)*tpxnum),&
                    &topoz(itpx+itpy*tpxnum),&
                    &topoz(itpx+1+itpy*tpxnum),&
                    &tpx1,tpy1,tpx2,tpy2,lontmp,lantmp)
        write(45,*)hrlen,elev
        write(50,*)hrlen,elev
    end do
end do
write(50,*)hrlen,0.000
write(50,*)0.000,0.000
close(45)
close(50)


stop
end



