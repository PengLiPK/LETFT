! This subroutine finds Frechet Derivative for inverse with 3d non-uniform
! rectangle velocity grid, in spherical coordinates, and earthquake 
! locations.
!-------------------------------------------------------------------------
subroutine frech_reg2sph3d2rl(vel,path,ipath,evnid,ttlevn,slevn,fd,layer,&
            &nl,topoxy,topoz,tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d) :: vel(maxgrid3d)
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: x,y,z,xlc,ylc,zlc
real(kind=8) :: length
real(kind=8) :: coflower(4),cofupper(4),zcof(2)
real(kind=8) :: tpminx,tpmaxx,tpminy,tpmaxy
real(kind=8) :: xtmp(4),ytmp(4),ztmp
real(kind=8) :: slevn,dlength,devnx,devny,devnz
integer :: lowerp(4),nlowerp
integer :: upperp(4),nupperp
integer :: nl(3)
integer :: n(8)
integer :: ipath
integer :: evnid,ttlevn
integer :: tpxnum,tpynum
integer :: i,j,k,updown


fd=0d0

! Calculate derivative of earthquake location.
call sphdist(path(1)%x,path(1)%y,path(1)%z,&
            &path(4)%x,path(4)%y,path(4)%z,&
            &dlength)
call sphdist(path(1)%x,path(1)%y,path(1)%z,&
            &path(4)%x,path(1)%y,path(1)%z,&
            &devnx)
call sphdist(path(1)%x,path(1)%y,path(1)%z,&
            &path(1)%x,path(4)%y,path(1)%z,&
            &devny)
devnz=path(4)%z-path(1)%z
if(path(4)%x .lt. path(1)%x)then
    devnx=-devnx
end if
if(path(4)%y .lt. path(1)%y)then
    devny=-devny
end if
fd(1+(evnid-1)*4)=1
fd(2+(evnid-1)*4)=-devnx*slevn/dlength
fd(3+(evnid-1)*4)=-devny*slevn/dlength
fd(4+(evnid-1)*4)=-devnz*slevn/dlength
!if(evnid .eq. 53)then
!!if(devnx .eq. 0d0 .and. devnx .eq. 0d0 .and. devnz .eq. 0d0)then
!    write(*,*)"error in strct.f90!",evnid
!    write(*,*)"dx,dy,dz,ds:",devnx,devny,devnz,dlength
!    write(*,*)slevn
!    write(*,*)path(1),path(4)
!    write(*,*)fd(209),fd(210),fd(211),fd(212)
!    stop
!end if

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0
    
    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)


    call locatcood3d2(n,nl(1),layer(1:maxgrd1d,1),nl(2),layer(1:maxgrd1d,2)&
                    &,nl(3),layer(1:maxgrd1d,3),x,y,z)
    
    ! Calculate the coefficiets of the nodes on the lower layer.
    ! Nodes in the air will not contribute.
    k=0
    do j=1,4
        call psurf(vel(n(j))%x,vel(n(j))%y,vel(n(j))%z,&
                  &topoxy,topoz,updown,tpminx,tpmaxx,&
                  &tpminy,tpmaxy,tpxnum,tpynum)
        if(updown .eq. 1)then
            k=k+1
            lowerp(k)=n(j)
        end if
    end do
    nlowerp=k

    ! Covert the earth coordinate to earth surface local Cartesian coordinate.
    ! Because *coef subroutines use Cartesian coordinate.
    if(nlowerp .gt. 1)then
        do j=1,nlowerp
            call ear2loc(xtmp(j),ytmp(j),ztmp,vel(lowerp(j))%x,&
                       &vel(lowerp(j))%y,vel(n(1))%z,&
                       &vel(n(1))%x,vel(n(1))%y,&
                       &vel(n(1))%z)
        end do
        call ear2loc(xlc,ylc,zlc,x,y,vel(n(1))%z,&
                    &vel(n(1))%x,vel(n(1))%y,vel(n(1))%z)
    end if

    coflower=0d0
    select case(nlowerp)
    case(4)
        call quadcoef(coflower,xlc,ylc,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
    case(3)
        call tricoef(coflower,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),&
                    &ytmp(2),xtmp(3),ytmp(3))
    case(2)
        call bicoef(coflower,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
    case(1)
        coflower(1)=1.0d0
    end select


    ! Calculate the coefficiets of the nodes on the upper layer.
    ! Nodes in the air will not contribute.
    k=0
    do j=5,8
        call psurf(vel(n(j))%x,vel(n(j))%y,vel(n(j))%z,&
                  &topoxy,topoz,updown,tpminx,tpmaxx,&
                  &tpminy,tpmaxy,tpxnum,tpynum)
        if(updown .eq. 1)then
            k=k+1
            upperp(k)=n(j)
        end if
    end do
    nupperp=k

    ! Covert the earth coordinate to earth surface local Cartesian coordinate.
    ! Because *coef subroutines use Cartesian coordinate.
    ! xtmp, ytmp, ztmp are used the second time here.
    if(nupperp .gt. 1)then
        do j=1,nupperp
            call ear2loc(xtmp(j),ytmp(j),ztmp,vel(upperp(j))%x,&
                       &vel(upperp(j))%y,vel(n(8))%z,&
                       &vel(n(1))%x,vel(n(1))%y,&
                       &vel(n(1))%z)
        end do
        call ear2loc(xlc,ylc,zlc,x,y,vel(n(8))%z,&
                    &vel(n(1))%x,vel(n(1))%y,vel(n(1))%z)
    end if


    cofupper=0d0
    select case(nupperp)
    case(4)
        call quadcoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
    case(3)
        call tricoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),&
                    &ytmp(2),xtmp(3),ytmp(3))
    case(2)
        call bicoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
    case(1)
        cofupper(1)=1.0d0
    end select


    ! z-direction coefficients
    if(nupperp .eq. 0 .and. nlowerp .ne. 0)then
        zcof(1)=1.0d0
        zcof(2)=0d0
    else if(nlowerp .eq. 0 .and. nupperp .ne. 0)then
        zcof(1)=0d0
        zcof(2)=1.0d0
    else if(nlowerp .eq. 0 .and. nupperp .eq. 0)then
        zcof(1)=0d0
        zcof(2)=0d0
    else
        zcof(1)=(vel(n(8))%z-z)/(vel(n(8))%z-vel(n(1))%z)
        zcof(2)=(z-vel(n(1))%z)/(vel(n(8))%z-vel(n(1))%z)
    end if

    
    fd(n(1)+ttlevn*4)=length*zcof(1)*coflower(1)+fd(n(1))
    fd(n(2)+ttlevn*4)=length*zcof(1)*coflower(2)+fd(n(2))
    fd(n(3)+ttlevn*4)=length*zcof(1)*coflower(3)+fd(n(3))
    fd(n(4)+ttlevn*4)=length*zcof(1)*coflower(4)+fd(n(4))

    fd(n(5)+ttlevn*4)=length*zcof(2)*cofupper(1)+fd(n(5))
    fd(n(6)+ttlevn*4)=length*zcof(2)*cofupper(2)+fd(n(6))
    fd(n(7)+ttlevn*4)=length*zcof(2)*cofupper(3)+fd(n(7))
    fd(n(8)+ttlevn*4)=length*zcof(2)*cofupper(4)+fd(n(8))

end do


return
end subroutine frech_reg2sph3d2rl
!-----------------------------------------------------------------------

