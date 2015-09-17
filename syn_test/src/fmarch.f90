module fmarch

implicit none

contains

!----------------------------------------------------------------------------
! This subroutine uses 1st order upwind method to update the travel times
! of narrow band grids.
!-----------------------------------------------------------------------------
subroutine march1(num,travelt,s,xnum,ttlgnum,dzz)

!
!--------------------------Paramters---------------------------------------
! travelt:
!       Travel time and coordinates of FMM nodes.
! v:
!       Velocities of FMM nodes.
! dyy
!       Grid spacing of FMM in km unit of y direction.
! temptx, tempty
!       Min (T(i-1,j), T(i+1,j)), Min (T(i,j-1), T(i,j+1)).
! a, b, c: 
!       The coeffients of quadratic finite difference form of Eikonal equation.
! num:
!       The number of new alive nodes, used for updating the travel times of 
!       new narrow band nodes.
! xnum:
!       The number of nodes in x directions.
! ttlgnum:
!       Total numbers of nodes in study area.
! gnum:
!       The number of neighbours of the new alive nodes.
! iudlr1, iudlr2:
!       Counters to represents the neighbours of the new alive nodes.
! xsol, ysol:
!       Determine the forms of a, b and c.
!--------------------------------------------------------------------------

use strct
implicit none
type(tstrct3d) :: travelt(maxgrid3d)
real(kind=8) :: s(maxgrid3d)
real(kind=8) :: dyy,dzz
real(kind=8) :: temptx,tempty
real(kind=8) :: a,b,c
integer :: num
integer :: xnum
integer :: ttlgnum
integer :: gnum
integer :: iudlr1,iudlr2
integer :: xsol,ysol

xsol=0
ysol=0
do iudlr1=-1,1,2
  do iudlr2=-1,1,2
  ! Update right, left, down and up grid
     gnum=num+(iudlr1-1)*iudlr2/2+&
     &    xnum*((iudlr1+1)*iudlr2/2)
     ! Check if these grids exist.
     if((gnum .gt. 0) .and. (gnum .le. ttlgnum))then
        if(((dble(iudlr2)*travelt(gnum)%x) .lt. &
         &  (dble(iudlr2)*travelt(num)%x)) .or. &
         & abs(travelt(gnum)%x-travelt(num)%x) .lt.&
         & 1d-4)then
            
            ! Make sure the  grid is not alive
            if(travelt(gnum)%stat .lt. 1)then

                ! Check upwind direction of x coodinate
                if((travelt(gnum+1)%x .gt. travelt(gnum)%x) .and.&
                & (travelt(gnum+1)%stat .eq. 1))then
                    if((travelt(gnum-1)%x .lt. travelt(gnum)%x) .and.&
                    & (travelt(gnum-1)%stat .eq. 1))then
                        temptx=min(travelt(gnum-1)%t,travelt(gnum+1)%t)
                        xsol=1
                    else
                        temptx=travelt(gnum+1)%t
                        xsol=1
                    end if
                else if((travelt(gnum-1)%x .lt. travelt(gnum)%x)&
                & .and. (travelt(gnum-1)%stat .eq. 1))then
                    temptx=travelt(gnum-1)%t
                    xsol=1
                else
                    temptx=0d0
                    xsol=0
                end if

                ! Check upwind direction of y coodinate
                if(((gnum-xnum) .gt. 0) .and. &
                & (travelt(gnum-xnum)%stat .eq. 1))then
                    if(((gnum+xnum) .le. ttlgnum) .and. &
                    & (travelt(gnum+xnum)%stat .eq. 1))then
                        tempty=min(travelt(gnum-xnum)%t,&
                        &      travelt(gnum+xnum)%t)
                        ysol=1
                    else
                        tempty=travelt(gnum-xnum)%t
                        ysol=1
                    end if
                else if(((gnum+xnum) .le. ttlgnum) .and. &
                &  (travelt(gnum+xnum)%stat .eq. 1))then
                    tempty=travelt(gnum+xnum)%t
                    ysol=1
                else
                    tempty=0d0
                    ysol=0
                end if

                dzz=1d0
                dyy=1d0

                a=dble(xsol)/(travelt(gnum)%dxx**2)+dble(ysol)/(dyy**2)
                b=-2.0d0*((dble(xsol)*temptx)/(travelt(gnum)%dxx**2)+&
                & (dble(ysol)*tempty)/(dyy**2))
                c=(dble(xsol)*((temptx/travelt(gnum)%dxx)**2))+ &
                & (dble(ysol)*((tempty/dyy)**2))-(s(gnum)**2)
                travelt(gnum)%t=(-b+sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
                travelt(gnum)%stat=2
            end if
        end if
     end if
  end do
end do

end subroutine march1
!---------------------------------------------------------------------------------
!---------------Subroutine march1 ends here---------------------------------------
!---------------------------------------------------------------------------------




!---------------------------------------------------------------------------------
! This subroutine uses 2nd order upwind method to update the travel times
! of narrow band grids.
!---------------------------------------------------------------------------------
subroutine march2(num,travelt,s,xnum,ynum,ttlgnum,dzz,nbnode,nbnum)

!
!--------------------------Paramters----------------------------------------------
! travelt:
!       Travel time and coordinates of FMM nodes.
! v:
!       Velocities of FMM nodes.
! temptx, tempty
!       Min (temptx1, temptx2), Min (tempty1, tempty2).
! temptx1, tempty1
!       2T(i-1,j)-0.5T(i-2,j); 2T(i,j-1)-0.5T(i,j-2).
! temptx2, tempty2
!       2T(i+1,j)-0.5T(i+2,j); 2T(i,j+1)-0.5T(i,j+2).
! a, b, c: 
!       The coeffients of quadratic finite difference form of Eikonal equation.
! sx, sy:
!       The weighted average velocity in x and y direction.
! savr:
!       The velocity used in finite difference form of Eikonal equation.
! nbnode:
!       The new narrow band node updated by subroutine "march2".
! nbnum:
!       The total number of new narrow band node.
! num:
!       The number of new alive nodes, used for updating the travel times of 
!       new narrow band nodes.
! xnum:
!       The number of nodes in x directions.
! ttlgnum:
!       Total numbers of nodes in study area.
! gnum:
!       The number of neighbours of the new alive nodes.
! iudlr1, iudlr2:
!       Counters to represents the neighbours of the new alive nodes.
! xsol, ysol, zsol:
!       Determine the forms of a, b and c.
! xsubsol, ysubsol, zsubsol:
!       Determine the forms of a, b and c.
! xmthd, ymthd, zmthd:
!       Determine the forms of a, b and c.
!--------------------------------------------------------------------------

use strct
implicit none
type(tstrct3d) :: travelt(maxgrid3d)
real(kind=8) :: s(maxgrid3d)
real(kind=8) :: dzz
real(kind=8) :: temptx,tempty,temptz
real(kind=8) :: temptx1,temptx2
real(kind=8) :: tempty1,tempty2
real(kind=8) :: temptz1,temptz2
real(kind=8) :: a,b,c
real(kind=8) :: bac
real(kind=8) :: sx,sy,sz
real(kind=8) :: savr
integer :: nbnode(maxnbnode3d)
integer :: nbnum
integer :: num
integer :: xnum,ynum
integer :: ttlgnum
integer :: gnum(6)
integer :: xsol,ysol,zsol
integer :: xsubsol,ysubsol,zsubsol
integer :: xmthd,ymthd,zmthd
integer :: ig,gmax
integer :: tmpgnum
integer :: is

nbnum=0
xsol=0
ysol=0
zsol=0


gmax=0
! Check neighbors in x direction are exist.
do ig=-1,1,2
    tmpgnum=num+ig
    if((tmpgnum .gt. 0) .and. (tmpgnum .le. ttlgnum) .and.&
      &(dble(ig)*(travelt(tmpgnum)%x-travelt(num)%x) .gt. 0d0))then
        gmax=gmax+1
        gnum(gmax)=tmpgnum
    end if
end do
    
! Check neighbors in y direction are exist or not.
do ig=-1,1,2
    tmpgnum=num+ig*xnum
    if((tmpgnum .gt. 0) .and. (tmpgnum .le. ttlgnum) .and.&
      &(dble(ig)*(travelt(tmpgnum)%y-travelt(num)%y) .gt. 0d0))then
        gmax=gmax+1
        gnum(gmax)=tmpgnum
    end if
end do
    
! Check neighbors in z direction are exist or not.
do ig=-1,1,2
    tmpgnum=num+ig*xnum*ynum
    if((tmpgnum .gt. 0) .and. (tmpgnum .le. ttlgnum))then
        gmax=gmax+1
        gnum(gmax)=tmpgnum
    end if
end do
    

do ig=1,gmax
    
    ! Make sure the  grid is not alive
    if(travelt(gnum(ig))%stat .lt. 1)then

        ! Check upwind direction of x coodinate
        if(((gnum(ig)+1) .gt. 0) .and. ((gnum(ig)+1) .le. ttlgnum .and.&
        & (travelt(gnum(ig)+1)%x .gt. travelt(gnum(ig))%x) .and.&
        & (travelt(gnum(ig)+1)%stat .eq. 1)))then
            if(((gnum(ig)-1) .gt. 0) .and. ((gnum(ig)-1) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-1)%x .lt. travelt(gnum(ig))%x) .and.&
            & (travelt(gnum(ig)-1)%stat .eq. 1)))then

                ! Check second order
                if(((gnum(ig)+2) .gt. 0) .and. ((gnum(ig)+2) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2)%x .gt. travelt(gnum(ig))%x) .and.&
                & (travelt(gnum(ig)+2)%stat .eq. 1)))then
                    if(((gnum(ig)-2) .gt. 0) .and. ((gnum(ig)-2) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2)%x .lt. travelt(gnum(ig))%x) .and.&
                    & (travelt(gnum(ig)-2)%stat .eq. 1)))then
                        xsol=2
                        xsubsol=-22
                    else
                    ! not x-2
                        xsol=2
                        xsubsol=-12
                    end if
                else
                ! not x+2
                    if(((gnum(ig)-2) .gt. 0) .and. ((gnum(ig)-2) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2)%x .lt. travelt(gnum(ig))%x) .and.&
                    & (travelt(gnum(ig)-2)%stat .eq. 1)))then
                        xsol=2
                        xsubsol=-21
                    else
                    ! not x-2
                        xsol=1
                        xsubsol=-11
                    end if
                end if
            else 
            ! not x-1
                if(((gnum(ig)+2) .gt. 0) .and. ((gnum(ig)+2) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2)%x .gt. travelt(gnum(ig))%x) .and.&
                & (travelt(gnum(ig)+2)%stat .eq. 1)))then
                    xsol=2
                    xsubsol=2
                else
                    xsol=1
                    xsubsol=1
                end if
            end if
        else 
        ! not x+1
            if(((gnum(ig)-1) .gt. 0) .and. ((gnum(ig)-1) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-1)%x .lt. travelt(gnum(ig))%x) .and.&
            & (travelt(gnum(ig)-1)%stat .eq. 1)))then
                if(((gnum(ig)-2) .gt. 0) .and. ((gnum(ig)-2) .le. ttlgnum .and.&
                & (travelt(gnum(ig)-2)%x .lt. travelt(gnum(ig))%x) .and.&
                & (travelt(gnum(ig)-2)%stat .eq. 1)))then
                    xsol=2
                    xsubsol=-2
                else
                    xsol=1
                    xsubsol=-1
                end if
            else
            ! not x-1
                temptx=0d0
                sx=s(gnum(ig))
                xsol=0
                xsubsol=0
            end if
        end if

        ! Check upwind direction of y coodinate.
        if(((gnum(ig)+xnum) .gt. 0) .and. ((gnum(ig)+xnum) .le. ttlgnum .and.&
        & (travelt(gnum(ig)+xnum)%y .gt. travelt(gnum(ig))%y) .and.&
        & (travelt(gnum(ig)+xnum)%stat .eq. 1)))then
            if(((gnum(ig)-xnum) .gt. 0) .and. ((gnum(ig)-xnum) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-xnum)%y .lt. travelt(gnum(ig))%y) .and.&
            & (travelt(gnum(ig)-xnum)%stat .eq. 1)))then

                !Check second order
                if(((gnum(ig)+2*xnum) .gt. 0) .and. ((gnum(ig)+2*xnum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2*xnum)%y .gt. travelt(gnum(ig))%y) .and.&
                & (travelt(gnum(ig)+2*xnum)%stat .eq. 1)))then
                    if(((gnum(ig)-2*xnum) .gt. 0) .and. ((gnum(ig)-2*xnum) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2*xnum)%y .lt. travelt(gnum(ig))%y) .and.&
                    & (travelt(gnum(ig)-2*xnum)%stat .eq. 1)))then
                        ysol=2
                        ysubsol=-22
                    else
                        ysol=2
                        ysubsol=-12
                    end if
                else
                ! not y+2
                    if(((gnum(ig)-2*xnum) .gt. 0) .and. ((gnum(ig)-2*xnum) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2*xnum)%y .lt. travelt(gnum(ig))%y) .and.&
                    & (travelt(gnum(ig)-2*xnum)%stat .eq. 1)))then
                        ysol=2
                        ysubsol=-21
                    else
                    ! not y-2
                        ysol=1
                        ysubsol=-11
                    end if
                end if
            else 
            ! not y-1
                if(((gnum(ig)+2*xnum) .gt. 0) .and. ((gnum(ig)+2*xnum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2*xnum)%y .gt. travelt(gnum(ig))%y) .and.&
                & (travelt(gnum(ig)+2*xnum)%stat .eq. 1)))then
                    ysol=2
                    ysubsol=2
                else
                    ysol=1
                    ysubsol=1
                end if
            end if
        else 
        ! not y+1
            if(((gnum(ig)-xnum) .gt. 0) .and. ((gnum(ig)-xnum) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-xnum)%y .lt. travelt(gnum(ig))%y) .and.&
            & (travelt(gnum(ig)-xnum)%stat .eq. 1)))then
                if(((gnum(ig)-2*xnum) .gt. 0) .and. ((gnum(ig)-2*xnum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)-2*xnum)%y .lt. travelt(gnum(ig))%y) .and.&
                & (travelt(gnum(ig)-2*xnum)%stat .eq. 1)))then
                    ysol=2
                    ysubsol=-2
                else
                    ysol=1
                    ysubsol=-1
                end if
            else
            ! not y-1
                tempty=0d0
                sy=s(gnum(ig))
                ysol=0
                ysubsol=0
            end if
        end if


        ! Check upwind direction of z coodinate.
        if(((gnum(ig)+xnum*ynum) .gt. 0) .and. ((gnum(ig)+xnum*ynum) .le. ttlgnum .and.&
        & (travelt(gnum(ig)+xnum*ynum)%stat .eq. 1)))then
            if(((gnum(ig)-xnum*ynum) .gt. 0) .and. ((gnum(ig)-xnum*ynum) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-xnum*ynum)%stat .eq. 1)))then

                !Check second order
                if(((gnum(ig)+2*xnum*ynum) .gt. 0) .and. ((gnum(ig)+2*xnum*ynum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2*xnum*ynum)%stat .eq. 1)))then
                    if(((gnum(ig)-2*xnum*ynum) .gt. 0) .and. ((gnum(ig)-2*xnum*ynum) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2*xnum*ynum)%stat .eq. 1)))then
                        zsol=2
                        zsubsol=-22
                    else
                    ! not z-2
                        zsol=2
                        zsubsol=-12
                    end if
                else
                ! not z+2
                    if(((gnum(ig)-2*xnum*ynum) .gt. 0) .and. ((gnum(ig)-2*xnum*ynum) .le. ttlgnum .and.&
                    & (travelt(gnum(ig)-2*xnum*ynum)%stat .eq. 1)))then
                        zsol=2
                        zsubsol=-21
                    else
                    ! not z-2
                        zsol=1
                        zsubsol=-11
                    end if
                end if
            else 
            ! not z-1
                if(((gnum(ig)+2*xnum*ynum) .gt. 0) .and. ((gnum(ig)+2*xnum*ynum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)+2*xnum*ynum)%stat .eq. 1)))then
                    zsol=2
                    zsubsol=2
                else
                    zsol=1
                    zsubsol=1
                end if
            end if
        else 
        ! not z+1
            if(((gnum(ig)-xnum*ynum) .gt. 0) .and. ((gnum(ig)-xnum*ynum) .le. ttlgnum .and.&
            & (travelt(gnum(ig)-xnum*ynum)%stat .eq. 1)))then
                if(((gnum(ig)-2*xnum*ynum) .gt. 0) .and. ((gnum(ig)-2*xnum*ynum) .le. ttlgnum .and.&
                & (travelt(gnum(ig)-2*xnum*ynum)%stat .eq. 1)))then
                    zsol=2
                    zsubsol=-2
                else
                ! not z-2
                    zsol=1
                    zsubsol=-1
                end if
            else
            ! not z-1
                temptz=0d0
                sz=s(gnum(ig))
                zsol=0
                zsubsol=0
            end if
        end if


        do is=1,2
            ! Now we know sols and subsols, the next step is using sols to
            ! determine finite difference method (1st order or 2nd order) 
            ! and temptx, tempty, temptz.
            if(xsol .le. 1)then
                ! 1st order method
                ! x direction
                if(xsubsol .eq. 0)then
                    temptx=0d0
                    sx=s(gnum(ig))
                    xmthd=0
                else if(xsubsol .eq. -11)then
                    temptx1=travelt(gnum(ig)-1)%t
                    temptx2=travelt(gnum(ig)+1)%t
                    if(temptx1 .le. temptx2)then
                        temptx=temptx1
                        sx=(s(gnum(ig))+s(gnum(ig)-1))/2.0d0
                    else
                        temptx=temptx2
                        sx=(s(gnum(ig))+s(gnum(ig)+1))/2.0d0
                    end if
                    xmthd=1
                else if(xsubsol .eq. -1)then
                    temptx=travelt(gnum(ig)-1)%t
                    sx=(s(gnum(ig))+s(gnum(ig)-1))/2.0d0
                    xmthd=1
                else if(xsubsol .eq. 1)then
                    temptx=travelt(gnum(ig)+1)%t
                    sx=(s(gnum(ig))+s(gnum(ig)+1))/2.0d0
                    xmthd=1
                else
                    write(*,*)"xsubsol is wrong!! xsubsol=",xsubsol
                end if
            else if(xsol .eq. 2)then
                ! 2nd or mixed order
                if(xsubsol .eq. -22)then
                    temptx1=(2.0d0*travelt(gnum(ig)-1)%t-&
                    & 0.5d0*travelt(gnum(ig)-2)%t)
                    temptx2=(2.0d0*travelt(gnum(ig)+1)%t-&
                    & 0.5d0*travelt(gnum(ig)+2)%t)
                    if(temptx1 .le. temptx2)then
                        temptx=temptx1
                        sx=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-1)+s(gnum(ig)-2))/8.0d0
                    else
                        temptx=temptx2
                        sx=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+1)+s(gnum(ig)+2))/8.0d0
                    end if
                    xmthd=2
                    xsubsol=-11
                else if((xsubsol .eq. -12) .or. (xsubsol .eq. -21))then
                    temptx1=travelt(gnum(ig)-1)%t
                    temptx2=travelt(gnum(ig)+1)%t
                    if(temptx1 .le. temptx2)then
                        temptx=temptx1
                        sx=(s(gnum(ig))+s(gnum(ig)-1))/2.0d0
                        xmthd=1
                    else
                        temptx=temptx2
                        sx=(s(gnum(ig))+s(gnum(ig)+1))/2.0d0
                        xmthd=1
                    end if
                    xsubsol=-11
                else if(xsubsol .eq. -2)then
                    temptx=(2.0d0*travelt(gnum(ig)-1)%t-&
                    & 0.5d0*travelt(gnum(ig)-2)%t)
                    sx=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-1)+s(gnum(ig)-2))/8.0d0
                    xmthd=2
                    xsubsol=-1
                else if(xsubsol .eq. 2)then
                    temptx=(2.0d0*travelt(gnum(ig)+1)%t-&
                    & 0.5d0*travelt(gnum(ig)+2)%t)
                    sx=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+1)+s(gnum(ig)+2))/8.0d0
                    xmthd=2
                    xsubsol=1
                else
                    write(*,*)"xsubsol is wrong!! xsubsol=",xsubsol
                end if
                xsol=1
            else
                write(*,*)"xsol is wrong!! xsol=",xsol
            end if
                
            if(ysol .le. 1)then
                ! y direction
                ! 1st order
                if(ysubsol .eq. 0)then
                    tempty=0d0
                    sy=s(gnum(ig))
                    ymthd=0
                else if(ysubsol .eq. -11)then
                    tempty1=travelt(gnum(ig)-xnum)%t
                    tempty2=travelt(gnum(ig)+xnum)%t
                    if(tempty1 .le. tempty2)then
                        tempty=tempty1
                        sy=(s(gnum(ig))+s(gnum(ig)-xnum))/2.0d0
                    else
                        tempty=tempty2
                        sy=(s(gnum(ig))+s(gnum(ig)+xnum))/2.0d0
                    end if
                    ymthd=1
                else if(ysubsol .eq. -1)then
                    tempty=travelt(gnum(ig)-xnum)%t
                    sy=(s(gnum(ig))+s(gnum(ig)-xnum))/2.0d0
                    ymthd=1
                else if(ysubsol .eq. 1)then
                    tempty=travelt(gnum(ig)+xnum)%t
                    sy=(s(gnum(ig))+s(gnum(ig)+xnum))/2.0d0
                    ymthd=1
                else
                    write(*,*)"ysubsol is wrong!! ysubsol=",ysubsol
                end if
            else if(ysol .eq. 2)then
                ! 2nd or mixed order
                if(ysubsol .eq. -22)then
                    tempty1=(2.0d0*travelt(gnum(ig)-xnum)%t-&
                    & 0.5d0*travelt(gnum(ig)-2*xnum)%t)
                    tempty2=(2.0d0*travelt(gnum(ig)+xnum)%t-&
                    & 0.5d0*travelt(gnum(ig)+2*xnum)%t)
                    if(tempty1 .le. tempty2)then
                        tempty=tempty1
                        sy=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-xnum)+&
                        & s(gnum(ig)-2*xnum))/8.0d0
                    else
                        tempty=tempty2
                        sy=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+xnum)+&
                        & s(gnum(ig)+2*xnum))/8.0d0
                    end if
                    ymthd=2
                    ysubsol=-11
                else if((ysubsol .eq. -12) .or. (ysubsol .eq. -21))then
                    tempty1=travelt(gnum(ig)-xnum)%t
                    tempty2=travelt(gnum(ig)+xnum)%t
                    if(tempty1 .le. tempty2)then
                        tempty=tempty1
                        sy=(s(gnum(ig))+s(gnum(ig)-xnum))/2.0d0
                        ymthd=1
                    else
                        tempty=tempty2
                        sy=(s(gnum(ig))+s(gnum(ig)+xnum))/2.0d0
                        ymthd=1
                    end if
                    ysubsol=-11
                else if(ysubsol .eq. -2)then
                    tempty=(2.0d0*travelt(gnum(ig)-xnum)%t-&
                    & 0.5d0*travelt(gnum(ig)-2*xnum)%t)
                    sy=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-xnum)+&
                    & s(gnum(ig)-2*xnum))/8.0d0
                    ymthd=2
                    ysubsol=-1
                else if(ysubsol .eq. 2)then
                    tempty=(2.0d0*travelt(gnum(ig)+xnum)%t-&
                    & 0.5d0*travelt(gnum(ig)+2*xnum)%t)
                    sy=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+xnum)+&
                    & s(gnum(ig)+2*xnum))/8.0d0
                    ymthd=2
                    ysubsol=1
                else
                    write(*,*)"ysubsol is wrong!! ysubsol=",ysubsol
                end if
                ysol=1
            else
                write(*,*)"ysol is wrong!! ysol=",ysol
            end if
                

            if(zsol .le. 1)then
                ! z direction
                ! 1st order
                if(zsubsol .eq. 0)then
                    temptz=0d0
                    sz=s(gnum(ig))
                    zmthd=0
                else if(zsubsol .eq. -11)then
                    temptz1=travelt(gnum(ig)-xnum*ynum)%t
                    temptz2=travelt(gnum(ig)+xnum*ynum)%t
                    if(temptz1 .le. temptz2)then
                        temptz=temptz1
                        sz=(s(gnum(ig))+s(gnum(ig)-xnum*ynum))/2.0d0
                    else
                        temptz=temptz2
                        sz=(s(gnum(ig))+s(gnum(ig)+xnum*ynum))/2.0d0
                    end if
                    zmthd=1
                else if(zsubsol .eq. -1)then
                    temptz=travelt(gnum(ig)-xnum*ynum)%t
                    sz=(s(gnum(ig))+s(gnum(ig)-xnum*ynum))/2.0d0
                    zmthd=1
                else if(zsubsol .eq. 1)then
                    temptz=travelt(gnum(ig)+xnum*ynum)%t
                    sz=(s(gnum(ig))+s(gnum(ig)+xnum*ynum))/2.0d0
                    zmthd=1
                else
                    write(*,*)"zsubsol is wrong!! zsubsol=",zsubsol
                end if
            else if(zsol .eq. 2)then
                ! 2nd or mixed order
                if(zsubsol .eq. -22)then
                    temptz1=(2.0d0*travelt(gnum(ig)-xnum*ynum)%t-&
                    & 0.5d0*travelt(gnum(ig)-2*xnum*ynum)%t)
                    temptz2=(2.0d0*travelt(gnum(ig)+xnum*ynum)%t-&
                    & 0.5d0*travelt(gnum(ig)+2*xnum*ynum)%t)
                    if(temptz1 .le. temptz2)then
                        temptz=temptz1
                        sz=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-xnum*ynum)+&
                        & s(gnum(ig)-2*xnum*ynum))/8.0d0
                    else
                        temptz=temptz2
                        sz=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+xnum*ynum)+&
                        & s(gnum(ig)+2*xnum*ynum))/8.0d0
                    end if
                    zmthd=2
                    zsubsol=-11
                else if((zsubsol .eq. -12) .or. (zsubsol .eq. -21))then
                    temptz1=travelt(gnum(ig)-xnum*ynum)%t
                    temptz2=travelt(gnum(ig)+xnum*ynum)%t
                    if(temptz1 .le. temptz2)then
                        temptz=temptz1
                        sz=(s(gnum(ig))+s(gnum(ig)-xnum*ynum))/2.0d0
                        zmthd=1
                    else
                        temptz=temptz2
                        sz=(s(gnum(ig))+s(gnum(ig)+xnum*ynum))/2.0d0
                        zmthd=1
                    end if
                    zsubsol=-11
                else if(zsubsol .eq. -2)then
                    temptz=(2.0d0*travelt(gnum(ig)-xnum*ynum)%t-&
                    & 0.5d0*travelt(gnum(ig)-2*xnum*ynum)%t)
                    sz=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)-xnum*ynum)+&
                    & s(gnum(ig)-2*xnum*ynum))/8.0d0
                    zmthd=2
                    zsubsol=-1
                else if(zsubsol .eq. 2)then
                    temptz=(2.0d0*travelt(gnum(ig)+xnum*ynum)%t-&
                    & 0.5d0*travelt(gnum(ig)+2*xnum*ynum)%t)
                    sz=(3.0d0*s(gnum(ig))+4.0d0*s(gnum(ig)+xnum*ynum)+&
                    & s(gnum(ig)+2*xnum*ynum))/8.0d0
                    zmthd=2
                    zsubsol=1
                else
                    write(*,*)"zsubsol is wrong!! zsubsol=",zsubsol
                end if
                zsol=1
            else
                write(*,*)"zsol is wrong!! zsol=",zsol
            end if



            ! Calculate travel time
            if(xmthd .gt. 0)then
                if(ymthd .gt. 0)then
                    if(zmthd .gt. 0)then
                        savr=(sx+sy+sz)/3.0d0
                    else
                        savr=(sx+sy)/2.0d0
                    end if
                else
                    if(zmthd .gt. 0)then
                        savr=(sx+sz)/2.0d0
                    else
                        savr=sx
                    end if
                end if
            else 
                if(ymthd .gt. 0)then
                    if(zmthd .gt. 0)then
                        savr=(sy+sz)/2.0d0
                    else
                        savr=sy
                    end if
                else
                    if(zmthd .gt. 0)then
                        savr=sz
                    else
                        write(*,*)"xmthd,ymthd,zmthd:",xmthd,ymthd,zmthd,"error in march2!"
                    end if
                end if
            end if

            !write(*,*)"line1743",gnum(ig),xmthd,ymthd,zmthd
            if(xmthd .lt. 2)then
                if(ymthd .lt. 2)then
                    if(zmthd .lt. 2)then
                        a=dble(xmthd)/(travelt(gnum(ig))%dxx**2)&
                        &+dble(ymthd)/(travelt(gnum(ig))%dyy**2)&
                        &+dble(zmthd)/(dzz**2)
                        b=-2.0d0*((dble(xmthd)*temptx)/(travelt(gnum(ig))%dxx**2)&
                        &+(dble(ymthd)*tempty)/(travelt(gnum(ig))%dyy**2)&
                        &+(dble(zmthd)*temptz)/(dzz**2))
                        c=dble(xmthd)*((temptx/travelt(gnum(ig))%dxx)**2)&
                        &+dble(ymthd)*((tempty/travelt(gnum(ig))%dyy)**2)&
                        &+dble(zmthd)*((temptz/dzz)**2)&
                        &-(savr**2)
                    else if(zmthd .eq. 2)then
                        a=dble(xmthd)/(travelt(gnum(ig))%dxx**2)&
                        &+dble(ymthd)/(travelt(gnum(ig))%dyy**2)&
                        &+9.0d0/(4.0d0*(dzz**2))
                        b=-2.0d0*((dble(xmthd)*temptx)/(travelt(gnum(ig))%dxx**2)&
                        &+(dble(ymthd)*tempty)/(travelt(gnum(ig))%dyy**2))&
                        &-3.0d0*temptz/(dzz**2)
                        c=dble(xmthd)*((temptx/travelt(gnum(ig))%dxx)**2)&
                        &+dble(ymthd)*((tempty/travelt(gnum(ig))%dyy)**2)&
                        &+(temptz/dzz)**2&
                        &-(savr**2)
                    else
                        write(*,*)"zmthd:",zmthd,"error in march2!"
                    end if
                else if(ymthd .eq. 2)then
                    if(zmthd .lt. 2)then
                        a=dble(xmthd)/(travelt(gnum(ig))%dxx**2)&
                        &+9.0d0/(4.0d0*(travelt(gnum(ig))%dyy**2))&
                        &+dble(zmthd)/(dzz**2)
                        b=-2.0d0*((dble(xmthd)*temptx)/(travelt(gnum(ig))%dxx**2)&
                        &+(dble(zmthd)*temptz)/(dzz**2))&
                        &-3.0d0*tempty/(travelt(gnum(ig))%dyy**2)
                        c=dble(xmthd)*((temptx/travelt(gnum(ig))%dxx)**2)&
                        &+(tempty/travelt(gnum(ig))%dyy)**2&
                        &+dble(zmthd)*((temptz/dzz)**2)&
                        &-(savr**2)
                     else if(zmthd .eq. 2)then
                        a=dble(xmthd)/(travelt(gnum(ig))%dxx**2)&
                        &+9.0d0/(4.0d0*(travelt(gnum(ig))%dyy**2))&
                        &+9.0d0/(4.0d0*(dzz**2))
                        b=-2.0d0*(dble(xmthd)*temptx)/(travelt(gnum(ig))%dxx**2)&
                        &-3.0d0*temptz/(dzz**2)&
                        &-3.0d0*tempty/(travelt(gnum(ig))%dyy**2)
                        c=dble(xmthd)*((temptx/travelt(gnum(ig))%dxx)**2)&
                        &+(tempty/travelt(gnum(ig))%dyy)**2&
                        &+(temptz/dzz)**2&
                        &-(savr**2)
                    else
                        write(*,*)"zmthd:",zmthd,"error in march2!"
                    end if
                else
                    write(*,*)"ymthd:",ymthd,"error in march2!"
                end if
            else if(xmthd .eq. 2)then 
                if(ymthd .lt. 2)then
                    if(zmthd .lt. 2)then
                        a=9.0d0/(4.0d0*(travelt(gnum(ig))%dxx**2))&
                        &+dble(ymthd)/(travelt(gnum(ig))%dyy**2)&
                        &+dble(zmthd)/(dzz**2)
                        b=-3.0d0*temptx/(travelt(gnum(ig))%dxx**2)&
                        &-2.0d0*(dble(zmthd)*temptz)/(dzz**2)&
                        &-2.0d0*(dble(ymthd)*tempty)/(travelt(gnum(ig))%dyy**2)
                        c=(temptx/travelt(gnum(ig))%dxx)**2&
                        &+dble(ymthd)*((tempty/travelt(gnum(ig))%dyy)**2)&
                        &+dble(zmthd)*((temptz/dzz)**2)&
                        &-(savr**2)
                      !  write(*,*)"line1794",gnum(ig),xmthd,ymthd,zmthd,a,b,c,temptx,tempty,temptz
                     else if(zmthd .eq. 2)then
                        a=9.0d0/(4.0d0*(travelt(gnum(ig))%dxx**2))&
                        &+dble(ymthd)/(travelt(gnum(ig))%dyy**2)&
                        &+9.0d0/(4.0d0*(dzz**2))
                        b=-3.0d0*temptx/(travelt(gnum(ig))%dxx**2)&
                        &-3.0d0*temptz/(dzz**2)&
                        &-2.0d0*(dble(ymthd)*tempty)/(travelt(gnum(ig))%dyy**2)
                        c=(temptx/travelt(gnum(ig))%dxx)**2&
                        &+dble(ymthd)*((tempty/travelt(gnum(ig))%dyy)**2)&
                        &+(temptz/dzz)**2&
                        &-(savr**2)
                    else
                        write(*,*)"zmthd:",zmthd,"error in march2!"
                    end if
                else if(ymthd .eq. 2)then
                    if(zmthd .lt. 2)then
                        a=9.0d0/(4.0d0*(travelt(gnum(ig))%dxx**2))&
                        &+9.0d0/(4.0d0*(travelt(gnum(ig))%dyy**2))&
                        &+dble(zmthd)/(dzz**2)
                        b=-3.0d0*temptx/(travelt(gnum(ig))%dxx**2)&
                        &-2.0d0*(dble(zmthd)*temptz)/(dzz**2)&
                        &-3.0d0*tempty/(travelt(gnum(ig))%dyy**2)
                        c=(temptx/travelt(gnum(ig))%dxx)**2&
                        &+(tempty/travelt(gnum(ig))%dyy)**2&
                        &+dble(zmthd)*((temptz/dzz)**2)&
                        &-(savr**2)
                     else if(zmthd .eq. 2)then
                        a=9.0d0/(4.0d0*(travelt(gnum(ig))%dxx**2))&
                        &+9.0d0/(4.0d0*(travelt(gnum(ig))%dyy**2))&
                        &+9.0d0/(4.0d0*(dzz**2))
                        b=-3.0d0*temptx/(travelt(gnum(ig))%dxx**2)&
                        &-3.0d0*temptz/(dzz**2)&
                        &-3.0d0*tempty/(travelt(gnum(ig))%dyy**2)
                        c=(temptx/travelt(gnum(ig))%dxx)**2&
                        &+(tempty/travelt(gnum(ig))%dyy)**2&
                        &+(temptz/dzz)**2&
                        &-(savr**2)
                    else
                        write(*,*)"zmthd:",zmthd,"error in march2!"
                    end if
                else
                    write(*,*)"ymthd:",ymthd,"error in march2!"
                end if
            else
                write(*,*)"ymthd:",ymthd,"error in march2!"
            end if
            bac=b**2-4.0d0*a*c

            ! If bac > 0, exit; or adopt 1st order and recalculate.
            if(bac .ge. 0d0)then
                exit
            end if
        end do

        if(is .eq. 2)then
            write(*,*)is
        end if

        if(bac .lt. 0d0)then
            write(*,*)"error in march2",is,bac,a,b,c
            write(*,*)"subsol",xsubsol,ysubsol,zsubsol
            write(*,*)"savr,sx,sy,sz",savr,sx,sy,sz
            write(*,*)travelt(gnum(ig))
            stop
        end if
        
        travelt(gnum(ig))%t=(-b+sqrt(b**2-4.0d0*a*c))/(2.0d0*a)
        travelt(gnum(ig))%stat=2
        nbnum=nbnum+1
        nbnode(nbnum)=gnum(ig)

        if(travelt(gnum(ig))%t .gt. 1d5)then
            write(*,*)"error in march2",bac,a,b,c,xmthd,ymthd,zmthd
            write(*,*)"subsol",xsubsol,ysubsol,zsubsol
            write(*,*)"savr,sx,sy,sz",savr,sx,sy,sz
            write(*,*)travelt(gnum(ig))
            stop
        end if
    end if
end do

return
end subroutine march2

!-----------------------------------------------------------------------------------------
!---------------Subrouting march2 ends here-----------------------------------------------
!-----------------------------------------------------------------------------------------

end module fmarch
