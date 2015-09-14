!function: ** linterp (linear interpolation)
!          ** bilinear (bilinear interpolation)
!          ** trisurfintp (trisurfintp interpolation)
!          ** deg2rad (Degrees --> radians)
!          ** rad2deg (Degrees <-- radians)
!          ** trilinear (Trilinear interpolation)
!          ** trilinear2 (Trilinear interpolation with considering 
!                         topography)
!          ** ear3dlinear (Trilinear interpolation for earth coordinate)
!subroutine: ** nb_del
!            ** nb_del3d
!            ** nb_output
!            ** nb_output3d
!            ** bil_gradient (gradient of bilinear interpolation)
!            ** tril_gradient (gradient of Trilinear interpolation)
!            ** ear3dl_gradient (gradient of 3D linear interpolation in
!                                earth coordinate)
!            ** baryc (barycentric interpolation in 2D (triangle))
!            ** baryc3dxyz (barycentric interpolation in 3D (tetrahedron),
!                           cartisian coordinates.)
!            ** baryc3dsph (barycentric interpolation in 3D (tetrahedron))
!                           spherical coordinates.)
!            ** vorointp3d (Voronoi diagram interpolation in 3D (tetrahedron))
!            ** baryc_grad (gradient of barycentric interpolation)
!            ** localcood (locate the center of line segment in
!                       regular grids)
!            ** localcood3d (locate the center of line segment in
!                      3D regular grids)
!            ** baryc_locat (locate the triangle for barycentric
!                      interpolation)
!            ** baryc_locat3d (locate the point (x,y,z) in which 
!                     tetrahedron)
!            ** raypath (find raypath when giving source and receiver)
!            ** raypath3d (find raypath in 3D when giving source and receiver)
!            ** quadcoef (calculate the coeficents of the interpolation of four
!                        rectangle nodes)
!            ** tricoef (calculate the coeficients of the interpolation of three
!                        triangle nodes)
!            ** bicoef (calculate the coeficients of the interpolation of two
!                       nodes)
!            ** frech_regular (calculate frechet derivative of regular grids)
!            ** frech_regsph3d (calculate frechet derivative of 3d regular grids,
!                              in spherical coordinates)
!            ** frech_regsph3d2 (calculate frechet derivative of 3d regular grids,
!                              in spherical coordinates, nodes in the air will
!                              not contribute.)
!            ** frech_reg2sph3d (calculate frechet derivative of non-uniform 
!                             regular grids, in spherical coordinates)
!            ** frech_reg2sph3d2 (calculate frechet derivative of non-uniform 
!                             regular grids, in spherical coordinates, nodes in
!                             the air will not contribute.)
!            ** frech_reg2sph3d2rl (calculate frechet derivative of non-uniform 
!                             regular grids, in spherical coordinates, nodes in
!                             the air will not contribute, derivative of
!                             earthquakes' location are included)
!            ** frech_tri (calculate frechet derivative of Delauney triangle grids)
!            ** frech_tethsph (calculate frechet derivative of 3D tetrahedro grids,
!                             in spherical coordinates)
!            ** frech_tethsph2 (calculate frechet derivative of 3D tetrahedro grids, the
!                           nodes in the air will not contribute,in spherical coordinates)
!            ** teth_vol (calculate volumn of tetrahedron in cartisian
!                         coordinates)
!            ** teth_volsph (calculate volumn of tetrahedron in spherical
!                           coordinates)
!            ** addnodes (add nodes in Delaunay triangles where ray 
!                    weights are high)
!            ** addnodes3d (add nodes in tetrahedros where ray 
!                    weights are high)
!            ** sphdist (calculate distance in spherical coordinate)
!            ** xyzdist (calculate distance in Cartesian coordinate)
!            ** sph2xyz (spherical coordinates --> cartisian coordinates)
!            ** xyz2sph (spherical coordinates <-- cartisian coordinates)
!            ** ear2xyz (earth coordinates --> cartisian coordinates)
!            ** xyz2ear (earth coordinates <-- cartisian coordinates)
!            ** ear2loc (earth coordinates --> earth surface local cartisian coordinates)
!            ** loc2ear (earth coordinates <-- earth surface local cartisian coordinates)
!            ** psurf (determine whethter a point is above or under a surface
!                   , +z is above)
!            ** surfele (determine the elevation of a point on a surface)
!            ** mergefile (merge a list of files into a single file)
!            ** mergebf (merge a list of binary files into a single file)
!            ** directionvector (calculate the direction vector of a line)
!            ** vfrom1d (calculate velocity of a point from 1d velocity model)
!            ** vfrom1dct (calculate velocity of a point from 1d continous velocity model)
!            ** det3 (calculate determinant of a 3x3 matrix, n!=6)
!            ** det4 (calculate determinant of a 4x4 matrix, n!=24)
!            ** dllocat3d (Delaunay search for tetrahran, locate the tetrahedron
!                   of a given point)
!            ** updateheap
!            ** upheap
!            ** downheap
!            ** detlayer (Determine r in which layer)
!            ** locatcood3d2 (Locate p(x,y,z) in the non-uniform cubics)
!            ** wmatrix (Determine the weighting matrix from a residual vector)
!--------------------------------------------------------------------------
module strct

    implicit none
! Paramters for 2D.
    integer,parameter :: maxgrid=100000
    integer,parameter :: maxpathnode=10000
    integer,parameter :: maxtri=5000
    integer,parameter :: maxsource=500
    integer,parameter :: maxreceiver=500
    integer,parameter :: maxdata=25000
    integer,parameter :: maxvel=1000
    integer,parameter :: maxnbnode=5000

! Paramters for 3D
    ! Number of nodes in fmm
    integer,parameter :: maxgrid3d=2000000
    ! Number of nodes on one ray path
    integer,parameter :: maxpathnode3d=10000
    ! Number of tetrahedrons from adaptive nodes
    integer,parameter :: maxteth=500000
    ! Number of sources of fmm
    integer,parameter :: maxsource3d=100
    ! Number of travel times
    integer,parameter :: maxdata3d=1200000
    ! Number of velocity parameters
    integer,parameter :: maxvel3d=50000
    ! Number of narrow band nodes in fmm
    integer,parameter :: maxnbnode3d=100000
    ! Number of total events
    integer,parameter :: maxevn3d=100000
    ! Number of nodes in 1 dimension
    integer,parameter :: maxgrd1d=1000
   
! Parameters for mpi
    integer,parameter :: maxproc=128

! Parameters for physicis constants.
! radii=earth radii + elevation of highest point.
    real(kind=8),parameter :: radii=6.375079d3
    real(kind=8),parameter :: pi=3.14159d0

! Source data structure, 2D and 3D
    type srstrct
        real(kind=8) :: x
        real(kind=8) :: y
    end type

    type srstrct3d
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: z
    end type

! Nodes structure, 2D and 3D
! stat=1,0,-1 means alive, narrowband, far away.
! stat=2 is temporary, to avoid duplication.
! nbstat=0 mean the nodes are not in narrowband.        |
! nbstat>0 mean the nodes are in narrowband, in this cas|
! nbstat is the index of the narrowband heap
    type tstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: t
        real(kind=8) :: dxx
        integer :: num
        integer :: stat
    end type
       
    type tstrct3d
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: z
        real(kind=8) :: t
        real(kind=8) :: dxx
        real(kind=8) :: dyy
        integer :: num
        integer :: nbstat
        integer :: stat
    end type
       
    type pstrct
        type(tstrct),pointer :: p
    end type

    type pstrct3d
        type(tstrct3d),pointer :: p
    end type

    type nb_linklist
        type(tstrct),pointer :: p
        type(nb_linklist),pointer :: prev
        type(nb_linklist),pointer :: next
    end type

    type nb_linklist3d
        type(tstrct3d),pointer :: p
        type(nb_linklist3d),pointer :: prev
        type(nb_linklist3d),pointer :: next
    end type

! Velocity structure, 2D and 3D
    type vstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: vel
    end type

    type vstrct3d
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: z
        real(kind=8) :: vel
    end type

! Receiver structure, 2D and 3D
    type rcstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: t
    end type

    type rcstrct3d
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: z
        real(kind=8) :: t
    end type

! G 1st norm structure
    type gnstrct
        real(kind=8) :: val
        integer :: num
    end type

    type pgnstrct
        type(gnstrct),pointer :: p
    end type
contains

! Delete element in nb_linklist.
!-------------------------------------------------------------
subroutine nb_del(item)
implicit none
type(nb_linklist),pointer :: item
type(nb_linklist),pointer :: prev,next

prev=>item%prev
next=>item%next
deallocate(item)
if(associated(prev))prev%next=>next
if(associated(next))next%prev=>prev
item=>next

return

end subroutine nb_del
!-------------------------------------------------------------


! Delete element in nb_linklist3d.
!-------------------------------------------------------------
subroutine nb_del3d(item)
implicit none
type(nb_linklist3d),pointer :: item
type(nb_linklist3d),pointer :: prev,next

prev=>item%prev
next=>item%next
deallocate(item)
if(associated(prev))prev%next=>next
if(associated(next))next%prev=>prev
item=>next

return

end subroutine nb_del3d
!-------------------------------------------------------------



! Output nb_linklist.
!--------------------------------------------------------------
subroutine nb_output(list,fname)
implicit none
type(nb_linklist),pointer :: list,p
character(len=70) :: fname

if(associated(list%prev) .eqv. .false.)then
    p=>list%next
else
    p=>list
end if

if(trim(fname) .eq. 'screen')then
    do while(associated(p))
        write(*,*)p%p
        p=>p%next
    end do
else
    open(900,file=fname,status='replace')
    do while(associated(p))
        write(900,*)p%p
        p=>p%next
    end do
end if


return

end subroutine nb_output
!---------------------------------------------------------------




! Output nb_linklist3d.
!--------------------------------------------------------------
subroutine nb_output3d(list,fname)
implicit none
type(nb_linklist3d),pointer :: list,p
character(len=70) :: fname

if(associated(list%prev) .eqv. .false.)then
    p=>list%next
else
    p=>list
end if

if(trim(fname) .eq. 'screen')then
    do while(associated(p))
        write(*,*)p%p
        p=>p%next
    end do
else
    open(900,file=fname,status='replace')
    do while(associated(p))
        write(900,*)p%p
        p=>p%next
    end do
end if


return

end subroutine nb_output3d
!---------------------------------------------------------------


! Linear interpolation for f(x)
!---------------------------------------------------------------
real(kind=8) function linterp(f1,f2,x1,x2,x)

implicit none
real(kind=8) :: f1,f2
real(kind=8) :: x1,x2
real(kind=8) :: x

if( (x2-x1) .le. 0d0 )then
    write(*,*)"error in function linear(strct.f90)"&
            &,(x2-x1)
end if

linterp=(f1*(x2-x) + f2*(x-x1))/(x2-x1)

return
end function linterp
!-----------------------------------------------------------------------



! Bilinear interpolation for f(x,y)
!---------------------------------------------------------------
real(kind=8) function bilinear(f11,f21,f12,f22,x1,y1,x2,y2,x,y)

implicit none
real(kind=8) :: f11,f21,f12,f22
real(kind=8) :: x1,x2,y1,y2
real(kind=8) :: x,y

if( ((x2-x1)*(y2-y1)) .le. 0d0 )then
    write(*,*)"error in function bilinear(strct.f90)"&
            &,((x2-x1)*(y2-y1))
end if

bilinear=(f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y)&
   &+f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1))&
   &/((x2-x1)*(y2-y1))

return
end function bilinear
!-----------------------------------------------------------------------


! Interpolation of a triangle, value is on the surface of tri.
!---------------------------------------------------------------
real(kind=8) function trisurfintp(f1,f2,f3,x1,y1,x2,y2,x3,y3,x,y)

implicit none
real(kind=8) :: f1,f2,f3
real(kind=8) :: x1,x2,x3,y1,y2,y3
real(kind=8) :: x,y,para

para=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

if( para .eq. 0d0 )then
    write(*,*)"error in function trisurfintp(strct.f90)"&
            &,para
    stop
end if

trisurfintp=((x-x1)*(y3-y1)*(f2-f1)&
            &+(y-y1)*(f3-f1)*(x2-x1)&
            &-(x-x1)*(y2-y1)*(f3-f1)&
            &-(y-y1)*(f2-f1)*(x3-x1))&
            &/para+f1

return
end function trisurfintp
!-----------------------------------------------------------------------


! Degrees to radian.
!------------------------------------------------------------------------
real(kind=8) function deg2rad(degree)

implicit none
real(kind=8) :: degree

deg2rad=degree*pi/1.8d2

return
end function deg2rad
!------------------------------------------------------------------------


! Degrees to radian.
!------------------------------------------------------------------------
real(kind=8) function rad2deg(degree)

implicit none
real(kind=8) :: degree

rad2deg=degree*1.8d2/pi

return
end function rad2deg
!------------------------------------------------------------------------

! Trilinear interpolation for f(x,y,z)
! crd:
!       Coordinate, 0 for Cartesian, 1 for spherical
!------------------------------------------------------------------------
real(kind=8) function trilinear(f111,f211,f121,f221,f112,f212,f122,f222,&
                            &x1,y1,z1,x2,y2,z2,x,y,z,crd)

implicit none
real(kind=8) :: f111,f211,f121,f221
real(kind=8) :: f112,f212,f122,f222
real(kind=8) :: x1,x2,y1,y2,z1,z2
real(kind=8) :: x,y,z
real(kind=8) :: xp,yp,zp
real(kind=8) :: xn,yn,zn
integer :: crd

if( ((x2-x1)*(y2-y1)*(z2-z1)) .le. 0d0 )then
    write(*,*)"error in function trilinear(strct.f90)",&
            &((x2-x1)*(y2-y1)*(z2-z1)),x1,x2,y1,y2,z1,z2,&
            &x,y,z
    stop
end if

if(crd .eq. 0)then
    xp=x2-x
    yp=y2-y
    zp=z2-z
    xn=x-x1
    yn=y-y1
    zn=z-z1
else if(crd .eq. 1)then
    xp=(x2-x)*cos(y1*pi/1.8d2)*pi*radii/1.8d2
    yp=(y2-y)*pi*radii/1.8d2
    zp=z2-z
    xn=(x-x1)*cos(y1*pi/1.8d2)*pi*radii/1.8d2
    yn=(y-y1)*pi*radii/1.8d2
    zn=z-z1
end if


trilinear=(((f111*xp + f211*xn)*yp&
          &+(f121*xp + f221*xn)*yn)*zp&
         &+((f112*xp + f212*xn)*yp&
          &+(f122*xp + f222*xn)*yn)*zn)&
          &/((xp+xn)*(yp+yn)*(zp+zn))

return
end function trilinear
!-------------------------------------------------------------------


! Trilinear interpolation for f(x,y,z) with considering topography.
! crd:
!       Coordinate, 0 for Cartesian, 1 for spherical
!------------------------------------------------------------------------
real(kind=8) function trilinear2(f111,f211,f121,f221,f112,f212,f122,f222,&
                            &x1,y1,z1,x2,y2,z2,x,y,z,topoxy,topoz,tpminx,&
                            &tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)

implicit none
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: tpminx,tpmaxx,tpminy,tpmaxy
real(kind=8) :: f111,f211,f121,f221
real(kind=8) :: f112,f212,f122,f222
real(kind=8) :: x1,x2,y1,y2,z1,z2
real(kind=8) :: x,y,z,xlc,ylc,zlc
real(kind=8) :: xnd(8),ynd(8),znd(8)
real(kind=8) :: coflower(4),cofupper(4),zcof(2)
real(kind=8) :: xtmp(4),ytmp(4),ztmp,cofsum
integer :: lowerp(4),nlowerp
integer :: upperp(4),nupperp
integer :: tpxnum,tpynum
integer :: j,k,updown


if( ((x2-x1)*(y2-y1)*(z2-z1)) .le. 0d0 )then
    write(*,*)"error in function trilinear(strct.f90)",&
            &((x2-x1)*(y2-y1)*(z2-z1)),x1,x2,y1,y2,z1,z2,&
            &x,y,z
    stop
end if

xnd(1)=x1
ynd(1)=y1
znd(1)=z1
xnd(2)=x2
ynd(2)=y1
znd(2)=z1
xnd(3)=x1
ynd(3)=y2
znd(3)=z1
xnd(4)=x2
ynd(4)=y2
znd(4)=z1

xnd(5)=x1
ynd(5)=y1
znd(5)=z2
xnd(6)=x2
ynd(6)=y1
znd(6)=z2
xnd(7)=x1
ynd(7)=y2
znd(7)=z2
xnd(8)=x2
ynd(8)=y2
znd(8)=z2


! Calculate the coefficiets of the nodes on the lower layer.
! Nodes in the air will not contribute.
k=0
do j=1,4
    call psurf(xnd(j),ynd(j),znd(j),&
              &topoxy,topoz,updown,tpminx,tpmaxx,&
              &tpminy,tpmaxy,tpxnum,tpynum)
    if(updown .eq. 1)then
        k=k+1
        lowerp(k)=j
    end if
end do
nlowerp=k

! Covert the earth coordinate to earth surface local Cartesian coordinate.
! Because *coef subroutines use Cartesian coordinate.
if(nlowerp .gt. 1)then
    do j=1,nlowerp
        call ear2loc(xtmp(j),ytmp(j),ztmp,xnd(lowerp(j)),&
                   &ynd(lowerp(j)),z1,&
                   &x1,y1,z1)
    end do
    ! Project the point to lower layer.
    call ear2loc(xlc,ylc,zlc,x,y,z1,x1,y1,z1)
end if


coflower=0d0
select case(nlowerp)
case(4)
    call quadcoef(coflower,xlc,ylc,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
case(3)
    call tricoef(coflower,xlc,ylc,xtmp(1),ytmp(1),&
                &xtmp(2),ytmp(2),xtmp(3),ytmp(3))
case(2)
    call bicoef(coflower,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
case(1)
    coflower(1)=1.0d0
end select

! Check the coefficient.
cofsum=0d0
!if(nlowerp .gt. 1)then
!    do j=1,nlowerp
!        cofsum=cofsum+coflower(j)
!        if(coflower(j) .gt. 1.001d0 .and. nlowerp .ne. 3)then
!            write(*,*)"Coefficient in trilinear2 is larger than 1!!"
!            write(*,*)"nlowerp,j,coflower:",nlowerp,j
!            write(*,*)(coflower(k),k=1,nlowerp)
!            write(*,*)xlc,ylc
!            write(*,*)xtmp(1),ytmp(1)
!            write(*,*)xtmp(2),ytmp(2)
!            write(*,*)xtmp(3),ytmp(3)
!            write(*,*)x,y,z
!            write(*,*)x1,y1,z1
!            write(*,*)x2,y2,z2
!            stop
!        end if
!    end do
!end if
!if(cofsum .lt. (1.0d0-1.0d-11) .and. cofsum .gt. 0d0)then
!    write(*,*)"cofsum of lowerp",cofsum
!    stop
!end if


! Calculate the coefficiets of the nodes on the upper layer.
! Nodes in the air will not contribute.
k=0
do j=5,8
    call psurf(xnd(j),ynd(j),znd(j),&
              &topoxy,topoz,updown,tpminx,tpmaxx,&
              &tpminy,tpmaxy,tpxnum,tpynum)
    if(updown .eq. 1)then
        k=k+1
        upperp(k)=j
    end if
end do
nupperp=k

! Covert the earth coordinate to earth surface local Cartesian coordinate.
! Because *coef subroutines use Cartesian coordinate.
if(nupperp .gt. 1)then
    do j=1,nupperp
        call ear2loc(xtmp(j),ytmp(j),ztmp,xnd(upperp(j)),&
                   &ynd(upperp(j)),z2,&
                   &x1,y1,z1)
    end do
    ! Project the point to lower layer.
    call ear2loc(xlc,ylc,zlc,x,y,z2,x1,y1,z1)
end if

cofupper=0d0
select case(nupperp)
case(4)
    call quadcoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
case(3)
    call tricoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),&
                &xtmp(2),ytmp(2),xtmp(3),ytmp(3))
case(2)
    call bicoef(cofupper,xlc,ylc,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
case(1)
    cofupper(1)=1.0d0
end select


! Check the coefficient.
!cofsum=0d0
!if(nupperp .gt. 1)then
!    do j=1,nupperp
!        cofsum=cofsum+cofupper(j)
!        if(cofupper(j) .gt. 1.001d0 .and. nlowerp .ne. 3)then
!            write(*,*)"Coefficient in trilinear2 is larger than 1!!"
!            write(*,*)"nupperp,j,cofupper:",nupperp,j
!            write(*,*)(cofupper(k),k=1,nupperp)
!            write(*,*)xlc,ylc
!            write(*,*)xtmp(1),ytmp(1)
!            write(*,*)xtmp(4),ytmp(4)
!            write(*,*)x,y,z
!            write(*,*)x1,y1,z1
!            write(*,*)x2,y2,z2
!            stop
!        end if
!    end do
!end if
!if(cofsum .lt. (1.0d0-1.0d-11) .and. cofsum .gt. 0d0)then
!    write(*,*)"cofsum of upperp",cofsum
!    stop
!end if


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
    zcof(1)=(z2-z)/(z2-z1)
    zcof(2)=(z-z1)/(z2-z1)
end if

trilinear2=(f111*coflower(1)+f211*coflower(2)&
          &+f121*coflower(3)+f221*coflower(4))*zcof(1)&
         &+(f112*cofupper(1)+f212*cofupper(2)&
          &+f122*cofupper(3)+f222*cofupper(4))*zcof(2)


return
end function trilinear2
!-------------------------------------------------------------------


! Trilinear interpolation for p(loni,lani,zi) in earth spherical coordinates
!------------------------------------------------------------------------
real(kind=8) function ear3dlinear(f111,f211,f121,f221,f112,f212,f122,f222,&
                            &lon1,lan1,z1,lon2,lan2,z2,loni,lani,zi)

implicit none
real(kind=8) :: f111,f211,f121,f221
real(kind=8) :: f112,f212,f122,f222
real(kind=8) :: lon1,lon2,lan1,lan2,z1,z2
real(kind=8) :: lani,loni,zi
real(kind=8) :: reflan,reflon,refz
real(kind=8) :: xp,yp,zp
real(kind=8) :: x(8),y(8),z(8)
real(kind=8) :: p12,p34,p1234
real(kind=8) :: p56,p78,p5678


if( ((lon2-lon1)*(lan2-lan1)*(z2-z1)) .le. 0d0 )then
        write(*,*)"error in function ear3dlinear(strct.f90)",&
                    &((lon2-lon1)*(lan2-lan1)*(z2-z1)),&
                    &lon1,lan1,z1,lon2,lan2,z2,loni,lani,zi
        stop
end if

! Calculate the coordinates of 8 points in local coordinate
reflon=(lon1+lon2)/2.0d0
reflan=(lan1+lan2)/2.0d0
refz=z1

call ear2loc(x(1),y(1),z(1),lon1,lan1,z1,reflon,reflan,refz)

x(2)=-x(1)
y(2)=y(1)
z(2)=z(1)


call ear2loc(x(3),y(3),z(3),lon1,lan2,z1,reflon,reflan,refz)

x(4)=-x(3)
y(4)=y(3)
z(4)=z(3)

call ear2loc(x(5),y(5),z(5),lon1,lan1,z2,reflon,reflan,refz)

x(6)=-x(5)
y(6)=y(5)
z(6)=z(5)

call ear2loc(x(7),y(7),z(7),lon1,lan2,z2,reflon,reflan,refz)

x(8)=-x(7)
y(8)=y(7)
z(8)=z(7)


! Covert the coordinates of interpolation point to local coordinates.

call ear2loc(xp,yp,zp,loni,lani,zi,reflon,reflan,refz)


p12=((x(2)-xp)*f111+(xp-x(1))*f211)/(x(2)-x(1))
p34=((x(4)-xp)*f121+(xp-x(3))*f221)/(x(4)-x(3))
p1234=((y(3)-yp)*p12+(yp-y(1))*p34)/(y(3)-y(1))

p56=((x(6)-xp)*f112+(xp-x(5))*f212)/(x(6)-x(5))
p78=((x(8)-xp)*f122+(xp-x(7))*f222)/(x(8)-x(7))
p5678=((y(7)-yp)*p12+(yp-y(5))*p34)/(y(7)-y(5))

ear3dlinear=((z(5)-zp)*p1234+(zp-z(1))*p5678)/(z(5)-z(1))

return
end function ear3dlinear
!-------------------------------------------------------------------



! Gradient of f(x,y) with bilinear interpolation
!-------------------------------------------------------------------
subroutine bil_gradient(gradx,grady,f11,f21,f12,f22,x1,y1,x2,y2,x,y)

implicit none
real(kind=8) :: f11,f21,f12,f22
real(kind=8) :: gradx,grady
real(kind=8) :: x1,x2,y1,y2
real(kind=8) :: x,y

if( ((x2-x1)*(y2-y1)) .le. 0d0 )then
    write(*,*)"error in subroutine bil_gradient(strct.f90)"&
            &,((x2-x1)*(y2-y1))
end if

gradx=((f21-f11)*(y2-y) + (f22-f12)*(y-y1))&
     &/((x2-x1)*(y2-y1))

grady=((f12-f11)*(x2-x) + (f22-f21)*(x-x1))&
     &/((x2-x1)*(y2-y1))

return

end subroutine bil_gradient
!----------------------------------------------------------------



! Gradient of f(x,y,z) with Trilinear interpolation
!-------------------------------------------------------------------
subroutine tril_gradient(gradx,grady,gradz,f111,f211,f121,f221,&
                &f112,f212,f122,f222,x1,y1,z1,x2,y2,z2,x,y,z)

implicit none
real(kind=8) :: f111,f211,f121,f221
real(kind=8) :: f112,f212,f122,f222
real(kind=8) :: x1,x2,y1,y2,z1,z2
real(kind=8) :: x,y,z
real(kind=8) :: xp,yp,zp
real(kind=8) :: xn,yn,zn
real(kind=8) :: prd
real(kind=8) :: gradx,grady,gradz


if( ((x2-x1)*(y2-y1)*(z2-z1)) .le. 0d0 )then
    write(*,*)"error in subroutine tril_gradient(strct.f90)"&
            &,((x2-x1)*(y2-y1)*(z2-z1))
    write(*,*)"x1,x2,y1,y2,z1,z2:"
    write(*,*)x1,x2,y1,y2,z1,z2
    stop
end if

xp=x2-x
yp=y2-y
zp=z2-z
xn=x-x1
yn=y-y1
zn=z-z1

prd=(x2-x1)*(y2-y1)*(z2-z1)

gradx=(((f211-f111)*yp + (f221-f121)*yn)*zp&
     &+((f212-f112)*yp + (f222-f122)*yn)*zn)&
     &/prd

grady=(((f121-f111)*xp + (f221-f211)*xn)*zp&
     &+((f122-f112)*xp + (f222-f212)*xn)*zn)&
     &/prd

gradz=(((f112-f111)*xp + (f212-f211)*xn)*yp&
     &+((f122-f121)*xp + (f222-f221)*xn)*yn)&
     &/prd


return

end subroutine tril_gradient
!----------------------------------------------------------------


! Trilinear interpolation for p(loni,lani,zi) in earth spherical coordinates
!------------------------------------------------------------------------
subroutine ear3dl_gradient(gx,gy,gz,f111,f211,f121,f221,f112,f212,f122,f222,&
                          &lon1,lan1,z1,lon2,lan2,z2,loni,lani,zi)

implicit none
real(kind=8) :: f111,f211,f121,f221
real(kind=8) :: f112,f212,f122,f222
real(kind=8) :: lon1,lon2,lan1,lan2,z1,z2
real(kind=8) :: lani,loni,zi
real(kind=8) :: reflan,reflon,refz
real(kind=8) :: xp,yp,zp
real(kind=8) :: x(8),y(8),z(8)
real(kind=8) :: p1,p2,p12
real(kind=8) :: p3,p4,p34
real(kind=8) :: gx,gy,gz

if( ((lon2-lon1)*(lan2-lan1)*(z2-z1)) .le. 0d0 )then
    write(*,*)"error in function ear3dl_gradient(strct.f90)",&
                &((lon2-lon1)*(lan2-lan1)*(z2-z1)),&
                &lon1,lan1,z1,lon2,lan2,z2,loni,lani,zi
    stop
end if

! Calculate the coordinates of 8 points in local coordinate
reflon=(lon1+lon2)/2.0d0
reflan=(lan1+lan2)/2.0d0
refz=z1

call ear2loc(x(1),y(1),z(1),lon1,lan1,z1,reflon,reflan,refz)

x(2)=-x(1)
y(2)=y(1)
z(2)=z(1)


call ear2loc(x(3),y(3),z(3),lon1,lan2,z1,reflon,reflan,refz)

x(4)=-x(3)
y(4)=y(3)
z(4)=z(3)

call ear2loc(x(5),y(5),z(5),lon1,lan1,z2,reflon,reflan,refz)

x(6)=-x(5)
y(6)=y(5)
z(6)=z(5)

call ear2loc(x(7),y(7),z(7),lon1,lan2,z2,reflon,reflan,refz)

x(8)=-x(7)
y(8)=y(7)
z(8)=z(7)


! Covert the coordinates of interpolation point to local coordinates.
call ear2loc(xp,yp,zp,loni,lani,zi,reflon,reflan,refz)


! Gradient in x direction
p1=(f211-f111)/(x(2)-x(1))
p2=(f221-f121)/(x(4)-x(3))
p12=((y(3)-yp)*p1+(yp-y(1))*p2)/(y(3)-y(1))

p3=(f212-f112)/(x(6)-x(5))
p4=(f222-f122)/(x(8)-x(7))
p34=((y(7)-yp)*p3+(yp-y(5))*p4)/(y(7)-y(5))

gx=((z(5)-zp)*p12+(zp-z(1))*p34)/(z(5)-z(1))


! Gradient in y direction
p1=((x(2)-xp)*f111+(xp-x(1))*f211)/(x(2)-x(1))
p2=((x(4)-xp)*f121+(xp-x(3))*f221)/(x(4)-x(3))
p12=(p2-p1)/(y(3)-y(1))

p3=((x(6)-xp)*f112+(xp-x(5))*f212)/(x(6)-x(5))
p4=((x(8)-xp)*f122+(xp-x(7))*f222)/(x(8)-x(7))
p34=(p4-p3)/(y(7)-y(5))

gy=((z(5)-zp)*p12+(zp-z(1))*p34)/(z(5)-z(1))


! Gradient in z direction
p12=((y(3)-yp)*p1+(yp-y(1))*p2)/(y(3)-y(1))

p34=((y(7)-yp)*p3+(yp-y(5))*p4)/(y(7)-y(5))

gz=(p34-p12)/(z(5)-z(1))

return
end subroutine ear3dl_gradient
!----------------------------------------------------------------------



! Barycentric interpolation in 2D (triangle)
!--------------------------------------------------------------------
subroutine baryc(f,f1,f2,f3,x1,y1,x2,y2,x3,y3,x,y)

implicit none
real(kind=8) :: f
real(kind=8) :: x,y
real(kind=8) :: f1,f2,f3
real(kind=8) :: x1,x2,x3
real(kind=8) :: y1,y2,y3
real(kind=8) :: area
real(kind=8) :: area1,area2,area3

area=abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

if( area .le. 0d0 )then
    write(*,*)"area of triangle is 0! error in subroutine &
            &baryc(strct.f90)",area
end if

area1=abs(x*y2-x2*y+x2*y3-x3*y2+x3*y-x*y3)
area2=abs(x1*y-x*y1+x*y3-x3*y+x3*y1-x1*y3)
area3=area-area1-area2

f=(f1*area1+f2*area2+f3*area3)/area

return

end subroutine baryc
!-------------------------------------------------------------------------



! Barycentric interpolation in tetrahedron with cartisian coordinates
!--------------------------------------------------------------------------
subroutine baryc3dxyz(f,f1,f2,f3,f4,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x,y,z)

implicit none
real(kind=8) :: f
real(kind=8) :: x,y,z
real(kind=8) :: f1,f2,f3,f4
real(kind=8) :: x1,x2,x3,x4
real(kind=8) :: y1,y2,y3,y4
real(kind=8) :: z1,z2,z3,z4
real(kind=8) :: vol
real(kind=8) :: vol1,vol2,vol3,vol4


call teth_vol(vol,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

if( vol .le. 0d0 )then
    write(*,*)"volume of triangle is 0! error in subroutine &
            &baryc3dxyz(strct.f90)",vol
end if

call teth_vol(vol1,x,y,z,x2,y2,z2,x3,y3,z3,x4,y4,z4)

call teth_vol(vol2,x1,y1,z1,x,y,z,x3,y3,z3,x4,y4,z4)

call teth_vol(vol3,x1,y1,z1,x2,y2,z2,x,y,z,x4,y4,z4)


vol4=vol-vol1-vol2-vol3


f=(f1*vol1+f2*vol2+f3*vol3+f4*vol4)/vol

return

end subroutine baryc3dxyz
!---------------------------------------------------------------------




! Barycentric interpolation in tetrahedron with cartisian coordinates
!--------------------------------------------------------------------------
subroutine baryc3dsph(f,f1,f2,f3,f4,lon1,lan1,z1,lon2,lan2,z2,lon3,lan3,z3,lon4,lan4,z4,lon,lan,z)

implicit none
real(kind=8) :: f
real(kind=8) :: lon,lan,z
real(kind=8) :: f1,f2,f3,f4
real(kind=8) :: lon1,lon2,lon3,lon4
real(kind=8) :: lan1,lan2,lan3,lan4
real(kind=8) :: z1,z2,z3,z4
real(kind=8) :: vol
real(kind=8) :: vol1,vol2,vol3,vol4


call teth_volsph(vol,lon1,lan1,z1,lon2,lan2,z2,lon3,lan3,z3,lon4,lan4,z4)

if( vol .le. 0d0 )then
    write(*,*)"volume of triangle is 0! error in subroutine &
            &baryc3dsph(strct.f90)",vol
end if

call teth_vol(vol1,lon,lan,z,lon2,lan2,z2,lon3,lan3,z3,lon4,lan4,z4)

call teth_vol(vol2,lon1,lan1,z1,lon,lan,z,lon3,lan3,z3,lon4,lan4,z4)

call teth_vol(vol3,lon1,lan1,z1,lon2,lan2,z2,lon,lan,z,lon4,lan4,z4)


vol4=vol-vol1-vol2-vol3


f=(f1*vol1+f2*vol2+f3*vol3+f4*vol4)/vol

return

end subroutine baryc3dsph
!---------------------------------------------------------------------


! Barycentric interpolation in 3D (tetrahedron)
!--------------------------------------------------------------------------
subroutine vorointp3dsph(f,f1,f2,f3,f4,&
                        &lon1,lan1,z1,lon2,lan2,z2,&
                        &lon3,lan3,z3,lon4,lan4,z4,&
                        &lon,lan,z)

implicit none
real(kind=8) :: f
real(kind=8) :: lon,lan,z
real(kind=8) :: f1,f2,f3,f4
real(kind=8) :: lon1,lon2,lon3,lon4
real(kind=8) :: lan1,lan2,lan3,lan4
real(kind=8) :: z1,z2,z3,z4
real(kind=8) :: sphd(4)
integer :: n


call sphdist(lon,lan,z,lon1,lan1,z1,sphd(1))
call sphdist(lon,lan,z,lon2,lan2,z2,sphd(2))
call sphdist(lon,lan,z,lon3,lan3,z3,sphd(3))
call sphdist(lon,lan,z,lon4,lan4,z4,sphd(4))
    
n=minloc(sphd,1)

select case (n)
    case (1)
        f=f1
    case (2)
        f=f2
    case (3)
        f=f3
    case (4)
        f=f4
end select


return

end subroutine vorointp3dsph
!---------------------------------------------------------------------


! Gradient of barycentric interpolation in 2D (triangle)
!--------------------------------------------------------------------
subroutine baryc_grad(gradx,grady,f1,f2,f3,x1,y1,x2,y2,x3,y3,x,y)

implicit none
real(kind=8) :: gradx,grady
real(kind=8) :: x,y
real(kind=8) :: f1,f2,f3
real(kind=8) :: x1,x2,x3
real(kind=8) :: y1,y2,y3
real(kind=8) :: area
real(kind=8) :: area1,area2,area3
real(kind=8) :: gx1,gx2,gx3
real(kind=8) :: gy1,gy2,gy3

area=abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

if( area .le. 0d0 )then
    write(*,*)"area of triangle is 0! error in subroutine &
            &baryc_grad(strct.f90)",area
end if

area1=x*y2-x2*y+x2*y3-x3*y2+x3*y-x*y3
if(area1 .ge. 0d0)then
    gx1=y2-x2*y+x2*y3-x3*y2+x3*y-y3
    gy1=x*y2-x2+x2*y3-x3*y2+x3-x*y3
else
    gx1=-(y2-x2*y+x2*y3-x3*y2+x3*y-y3)
    gy1=-(x*y2-x2+x2*y3-x3*y2+x3-x*y3)
end if

area2=x1*y-x*y1+x*y3-x3*y+x3*y1-x1*y3
if(area2 .ge. 0d0)then
    gx2=x1*y-y1+y3-x3*y+x3*y1-x1*y3
    gy2=x1-x*y1+x*y3-x3+x3*y1-x1*y3
else
    gx2=-(x1*y-y1+y3-x3*y+x3*y1-x1*y3)
    gy2=-(x1-x*y1+x*y3-x3+x3*y1-x1*y3)
end if

area3=x1*y2-x2*y1+x2*y-x*y2+x*y1-x1*y
if(area3 .ge. 0d0)then
    gx3=x1*y2-x2*y1+x2*y-y2+y1-x1*y
    gy3=x1*y2-x2*y1+x2-x*y2+x*y1-x1
else
    gx3=-(x1*y2-x2*y1+x2*y-y2+y1-x1*y)
    gy3=-(x1*y2-x2*y1+x2-x*y2+x*y1-x1)
end if

gradx=(f1*gx1+f2*gx2+f3*gx3)/area
grady=(f1*gy1+f2*gy2+f3*gy3)/area

return

end subroutine baryc_grad
!---------------------------------------------------------------------



! This subroutine determine the local coordinate of point (x,y).
!----------------------------------------------------------------------
subroutine localcood(x,y,travelt,prv,dx,dy,&
            &minx,maxx,miny,maxy,xnum,ynum)

implicit none
type(tstrct), target :: travelt(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: x,y
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: dx,dy
integer :: pnum
integer :: xnum,ynum
integer :: tmplnumx,tmplnumy

!write(*,*)"line574"
!write(*,*)size(prv)
!stop
!allocate(prv(1:4))
if(x .le. minx)then
    if(y .le. miny)then
        pnum=1
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+1
    else
        pnum=(ynum-2)*xnum+1
    end if
else if(x .lt. maxx)then
    if(y .le. miny)then
        tmplnumx=int((x-minx)/dx)
        pnum=tmplnumx+1
    else if(y .lt. maxy)then
        tmplnumx=int((x-minx)/dx)
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+tmplnumx+1
    else
        tmplnumx=int((x-minx)/dx)
        pnum=(ynum-2)*xnum+tmplnumx+1
    end if
else
    if(y .le. miny)then
        pnum=xnum-1
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+xnum-1
    else
        pnum=(ynum-2)*xnum+xnum-1
    end if
end if
prv(1)%p=>travelt(pnum)
prv(2)%p=>travelt(pnum+1)
prv(3)%p=>travelt(pnum+xnum)
prv(4)%p=>travelt(pnum+xnum+1)

!deallocate(prv)

return
end subroutine localcood
!----------------------------------------------------------------------



! This subroutine determine the locating grid of point (x,y,z).
!----------------------------------------------------------------------
subroutine localcood3d(x,y,z,travelt,prv,dx,dy,dz,&
            &minx,maxx,miny,maxy,minz,maxz,xnum,ynum,znum)

implicit none
type(tstrct3d), target :: travelt(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
real(kind=8) :: x,y,z
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: dx,dy,dz
integer :: pnum
integer :: xnum,ynum,znum
integer :: tmplnumx,tmplnumy,tmplnumz

!allocate(prv(1:8))

if(x .le. minx)then
    if(y .le. miny)then
        if(z .le. minz)then
            pnum=1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+1
        else
            pnum=(znum-2)*xnum*ynum+1
        end if
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        if(z .le. minz)then
            pnum=tmplnumy*xnum+1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+tmplnumy*xnum+1
        else
            pnum=(znum-2)*xnum*ynum+tmplnumy*xnum+1
        end if
    else
        if(z .le. minz)then
            pnum=(ynum-2)*xnum+1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+(ynum-2)*xnum+1
        else
            pnum=(znum-2)*xnum*ynum+(ynum-2)*xnum+1
        end if
    end if
else if(x .lt. maxx)then
    tmplnumx=int((x-minx)/dx)
    if(y .le. miny)then
        if(z .le. minz)then
            pnum=tmplnumx+1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+tmplnumx+1
        else
            pnum=(znum-2)*xnum*ynum+tmplnumx+1
        end if
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        if(z .le. minz)then
            pnum=tmplnumy*xnum+tmplnumx+1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+tmplnumy*xnum+tmplnumx+1
        else
            pnum=(znum-2)*xnum*ynum+tmplnumy*xnum+tmplnumx+1
        end if
    else
        if(z .le. minz)then
            pnum=(ynum-2)*xnum+tmplnumx+1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+(ynum-2)*xnum+tmplnumx+1
        else
            pnum=(znum-2)*xnum*ynum+(ynum-2)*xnum+tmplnumx+1
        end if
    end if
else
    if(y .le. miny)then
        if(z .le. minz)then
            pnum=xnum-1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+xnum-1
        else
            pnum=(znum-2)*xnum*ynum+xnum-1
        end if
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        if(z .le. minz)then
            pnum=tmplnumy*xnum+xnum-1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+tmplnumy*xnum+xnum-1
        else
            pnum=(znum-2)*xnum*ynum+tmplnumy*xnum+xnum-1
        end if
    else
        if(z .le. minz)then
            pnum=(ynum-2)*xnum+xnum-1
        else if(z .lt. maxz)then
            tmplnumz=int((z-minz)/dz)
            pnum=tmplnumz*xnum*ynum+(ynum-2)*xnum+xnum-1
        else
            pnum=(znum-2)*xnum*ynum+(ynum-2)*xnum+xnum-1
        end if
    end if
end if
prv(1)%p=>travelt(pnum)
prv(2)%p=>travelt(pnum+1)
prv(3)%p=>travelt(pnum+xnum)
prv(4)%p=>travelt(pnum+xnum+1)
prv(5)%p=>travelt(pnum+xnum*ynum)
prv(6)%p=>travelt(pnum+xnum*ynum+1)
prv(7)%p=>travelt(pnum+xnum*ynum+xnum)
prv(8)%p=>travelt(pnum+xnum*ynum+xnum+1)

!deallocate(prv)

return
end subroutine localcood3d
!----------------------------------------------------------------------



! This subroutine locates the triangle for barycentric interpolation
!----------------------------------------------------------------------
subroutine baryc_locat(x,y,velnode,tri,trinum,ptri)

implicit none
type(tstrct), target :: velnode(maxgrid)
type(pstrct), pointer :: ptri(:)
real(kind=8) :: x,y
real(kind=8) :: xt(3),yt(3)
real(kind=8) :: a,b,c
integer :: tri(maxtri,3)
integer :: trinum
integer :: it,i

!allocate(ptri(1:3))
do i=1,3
    nullify(ptri(i)%p)
end do


do it=1,trinum
    do i=1,3
        xt(i)=velnode(tri(it,i)+1)%x-x
        yt(i)=velnode(tri(it,i)+1)%y-y
    end do
    a=xt(1)*yt(2)-xt(2)*yt(1)
    b=xt(2)*yt(3)-xt(3)*yt(2)
    if( a*b .ge. -1D-10 )then
        c=xt(3)*yt(1)-xt(1)*yt(3)
        if( (b*c .ge. -1D-10) .and. (c*a .ge. -1D-10) )then
            ptri(1)%p=>velnode(tri(it,1)+1)
            ptri(2)%p=>velnode(tri(it,2)+1)
            ptri(3)%p=>velnode(tri(it,3)+1)
            exit
        end if
    end if
end do

if(associated(ptri(1)%p) .eqv. .false.)then
    write(*,*)"point",x,y,"is not in triangles"
    stop
end if

!deallocate(ptri)

return
end subroutine baryc_locat
!----------------------------------------------------------------------




! This subroutine locates the point (x,y,z) in tethahedrons.
!----------------------------------------------------------------------
subroutine baryc_locat3d(x,y,z,velnode,teth,tethnum,pteth,itet)

implicit none
type(tstrct3d), target :: velnode(maxgrid3d)
type(pstrct3d), pointer :: pteth(:)
real(kind=8) :: x,y,z
real(kind=8) :: x1,x2,x3,xp
real(kind=8) :: y1,y2,y3,yp
real(kind=8) :: z1,z2,z3,zp
real(kind=8) :: vol,vol1,vol2,vol3,vol4
integer :: teth(maxteth,4)
integer :: tethnum
integer :: it,i
integer :: itet

allocate(pteth(1:4))
do i=1,4
    nullify(pteth(i)%p)
end do

do it=1,tethnum
    
    x1=(velnode(teth(it,1)+1)%x-velnode(teth(it,4)+1)%x)*1.0d2
    y1=(velnode(teth(it,1)+1)%y-velnode(teth(it,4)+1)%y)*1.1d2
    z1=velnode(teth(it,1)+1)%z-velnode(teth(it,4)+1)%z
    x2=(velnode(teth(it,2)+1)%x-velnode(teth(it,4)+1)%x)*1.0d2
    y2=(velnode(teth(it,2)+1)%y-velnode(teth(it,4)+1)%y)*1.1d2
    z2=velnode(teth(it,2)+1)%z-velnode(teth(it,4)+1)%z
    x3=(velnode(teth(it,3)+1)%x-velnode(teth(it,4)+1)%x)*1.0d2
    y3=(velnode(teth(it,3)+1)%y-velnode(teth(it,4)+1)%y)*1.1d2
    z3=velnode(teth(it,3)+1)%z-velnode(teth(it,4)+1)%z
    
    call det3(vol,x1,y1,z1,x2,y2,z2,x3,y3,z3)
    
    xp=(x-velnode(teth(it,4)+1)%x)*1.0d2
    yp=(y-velnode(teth(it,4)+1)%y)*1.1d2
    zp=z-velnode(teth(it,4)+1)%z

    call det3(vol1,xp,yp,zp,x2,y2,z2,x3,y3,z3)


    if( (vol*vol1 .ge. 0d0) .or. (abs(vol1/vol) .lt. 1d-12) )then

        call det3(vol2,x1,y1,z1,xp,yp,zp,x3,y3,z3)

        if( (vol*vol2 .ge. 0d0) .or. (abs(vol2/vol) .lt. 1d-12) )then

            call det3(vol3,x1,y1,z1,x2,y2,z2,xp,yp,zp)

            if( (vol*vol3 .ge. 0d0) .or. (abs(vol3/vol) .lt. 1d-12) )then

                x1=(velnode(teth(it,1)+1)%x-x)*1.0d2
                y1=(velnode(teth(it,1)+1)%y-y)*1.1d2
                z1=velnode(teth(it,1)+1)%z-z
                x2=(velnode(teth(it,2)+1)%x-x)*1.0d2
                y2=(velnode(teth(it,2)+1)%y-y)*1.1d2
                z2=velnode(teth(it,2)+1)%z-z
                x3=(velnode(teth(it,3)+1)%x-x)*1.0d2
                y3=(velnode(teth(it,3)+1)%y-y)*1.1d2
                z3=velnode(teth(it,3)+1)%z-z
                
                call det3(vol4,x1,y1,z1,x2,y2,z2,x3,y3,z3)

                if( (vol*vol4 .ge. 0d0) .or. (abs(vol4/vol) .lt. 1d-12) )then
                    pteth(1)%p=>velnode(teth(it,1)+1)
                    pteth(2)%p=>velnode(teth(it,2)+1)
                    pteth(3)%p=>velnode(teth(it,3)+1)
                    pteth(4)%p=>velnode(teth(it,4)+1)
                    itet=it
                    exit
                end if
            end if
        end if
    end if
end do


if(associated(pteth(1)%p) .eqv. .false.)then
    write(*,*)"point",x,y,z,"is not in tethahedron"
    stop
end if

!deallocate(pteth)

return
end subroutine baryc_locat3d
!----------------------------------------------------------------------



! This subroutine finds ray path with knowing travel time of each grid.
!-----------------------------------------------------------------------
subroutine raypath(travelt,receiver,source,path,&
           &    ipth,dx,dy,minx,maxx,miny,maxy,xnum,ynum)

implicit none
type(srstrct) :: path(maxpathnode)
type(rcstrct) :: receiver
type(srstrct) :: source
type(tstrct), target :: travelt(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: dtx,dty
real(kind=8) :: dx,dy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
integer :: ipth
integer :: ith
integer :: xnum,ynum

allocate(prv(1:4))
ipth=1
path(ipth)%x=receiver%x
path(ipth)%y=receiver%y

do while((abs(path(ipth)%x-source%x) .gt. dx) .or.&
&       (abs(path(ipth)%y-source%y) .gt. dy))
        
! Find the local rectangle for next ray path point caculation.
!-------------------------------------------------------------------
    call localcood(path(ipth)%x,path(ipth)%y,travelt,prv,dx,dy,&
    &minx,maxx,miny,maxy,xnum,ynum)
!-------------------------------------------------------------------
! Local rectangle finding ends.

! Caculate next point of ray path
!-------------------------------------------------------------------

    call bil_gradient(dtx,dty,prv(1)%p%t,prv(2)%p%t,prv(3)%p%t,&
    &prv(4)%p%t,prv(1)%p%x,prv(1)%p%y,prv(4)%p%x,prv(4)%p%y,&
    &path(ipth)%x,path(ipth)%y)

    path(ipth+1)%x=path(ipth)%x-(((dx+dy)/8.0d0)*dtx/sqrt(dtx**2+dty**2))
    path(ipth+1)%y=path(ipth)%y-(((dx+dy)/8.0d0)*dty/sqrt(dtx**2+dty**2))
    
    ith=0
    if(path(ipth+1)%x .lt. minx)then
        path(ipth+1)%x=minx
        ith=ith+1
    else if(path(ipth+1)%x .gt. maxx)then
        path(ipth+1)%x=maxx
        ith=ith+1
    else if(path(ipth+1)%y .lt. miny)then
        path(ipth+1)%y=miny
        ith=ith+1
    else if(path(ipth+1)%y .gt. maxy)then
        path(ipth+1)%y=maxy
        ith=ith+1
    end if


    ipth=ipth+1
!------------------------------------------------------------------

end do

path(ipth)%x=source%x
path(ipth)%y=source%y
deallocate(prv)

return

end subroutine raypath
!----------------------------------------------------------------------




! This subroutine finds ray path with knowing travel time of each grid
! in 3D.
!-----------------------------------------------------------------------
subroutine raypath3d(travelt,receiver,source,path,ipth,dx,dy,dz,&
                &minx,maxx,miny,maxy,minz,maxz,xnum,ynum,znum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(rcstrct3d) :: receiver
type(srstrct3d) :: source
type(tstrct3d), target :: travelt(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
real(kind=8) :: dtx,dty,dtz
real(kind=8) :: dx,dy,dz
real(kind=8) :: x(2),y(2),z(2)
real(kind=8) :: pthx,pthy,pthz
real(kind=8) :: psegx,psegy,psegz
real(kind=8) :: nextpthx,nextpthy,nextpthz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: gsqrt
real(kind=8) :: reflen
integer :: ipth
integer :: ith
integer :: xnum,ynum,znum

allocate(prv(1:8))

ipth=1
path(ipth)%x=receiver%x
path(ipth)%y=receiver%y
path(ipth)%z=receiver%z

do while((abs(path(ipth)%x-source%x) .gt. dx) .or.&
&        (abs(path(ipth)%y-source%y) .gt. dy) .or.&
&        (abs(path(ipth)%z-source%z) .gt. dz))
        
! Find the local rectangle for next ray path point caculation.
!-------------------------------------------------------------------
    call localcood3d(path(ipth)%x,path(ipth)%y,path(ipth)%z,&
        &travelt,prv,dx,dy,dz,minx,maxx,miny,maxy,minz,maxz,&
        &xnum,ynum,znum)
!-------------------------------------------------------------------
! Local rectangle finding ends.

! Caculate next point of ray path
!-------------------------------------------------------------------

    x(1)=0d0
    y(1)=0d0
    z(1)=prv(1)%p%z
    
    call ear2loc(x(2),y(2),z(2),prv(8)%p%x,prv(8)%p%y,prv(8)%p%z,&
                &prv(1)%p%x,prv(1)%p%y,0d0)


    call ear2loc(pthx,pthy,pthz,path(ipth)%x,path(ipth)%y,path(ipth)%z,&
                &prv(1)%p%x,prv(1)%p%y,0d0)

    call tril_gradient(dtx,dty,dtz,prv(1)%p%t,prv(2)%p%t,prv(3)%p%t,&
        &prv(4)%p%t,prv(5)%p%t,prv(6)%p%t,prv(7)%p%t,prv(8)%p%t,&
        &x(1),y(1),z(1),x(2),y(2),z(2),&
        &pthx,pthy,pthz)


    gsqrt=sqrt(dtx**2+dty**2+dtz**2)

    reflen=(abs(x(2)-x(1))+abs(y(2)-y(1))+abs(z(2)-z(1)))/8.0d0
    psegx=reflen*dtx/gsqrt
    psegy=reflen*dty/gsqrt
    psegz=reflen*dtz/gsqrt

    pthx=pthx-psegx
    pthy=pthy-psegy
    pthz=pthz-psegz

    call loc2ear(nextpthx,nextpthy,nextpthz,&
                &prv(1)%p%x,prv(1)%p%y,0d0,pthx,pthy,pthz)

    path(ipth+1)%x=nextpthx
    path(ipth+1)%y=nextpthy
    path(ipth+1)%z=nextpthz

    ith=0
    if(path(ipth+1)%x .lt. minx)then
        path(ipth+1)%x=minx
        ith=ith+1
    else if(path(ipth+1)%x .gt. maxx)then
        path(ipth+1)%x=maxx
        ith=ith+1
    end if
    if(path(ipth+1)%y .lt. miny)then
        path(ipth+1)%y=miny
        ith=ith+1
    else if(path(ipth+1)%y .gt. maxy)then
        path(ipth+1)%y=maxy
        ith=ith+1
    end if
    if(path(ipth+1)%z .lt. minz)then
        path(ipth+1)%z=minz
        ith=ith+1
    else if(path(ipth+1)%z .gt. maxz)then
        path(ipth+1)%z=maxz
        ith=ith+1
    end if


    ipth=ipth+1
!------------------------------------------------------------------

end do

path(ipth)%x=source%x
path(ipth)%y=source%y
path(ipth)%z=source%z
deallocate(prv)

return

end subroutine raypath3d
!----------------------------------------------------------------------


! Calculate the coeficents of the interpolation of four nodes in a rectangle.
!-------------------------------------------------------------------------
subroutine quadcoef(coef,x,y,x1,y1,x2,y2)

implicit none
real(kind=8) :: coef(4)
real(kind=8) :: x1,y1,x2,y2
real(kind=8) :: x,y,a


a=(x2-x1)*(y2-y1)

coef(1)=(x2-x)*(y2-y)/a
coef(2)=(x-x1)*(y2-y)/a
coef(3)=(x2-x)*(y-y1)/a
coef(4)=(x-x1)*(y-y1)/a

return
end subroutine quadcoef
!-------------------------------------------------------------------------


! Calculate the coeficents of the interpolation or extapolation from 
! a triangle. Using plane method.
!-------------------------------------------------------------------------
subroutine tricoef(coef,x,y,x1,y1,x2,y2,x3,y3)

implicit none
real(kind=8) :: coef(4)
real(kind=8) :: x,y,x1,y1,x2,y2,x3,y3
real(kind=8) :: x01,x02,x03,y01,y02,y03
real(kind=8) :: tmpcoef


x01=x-x1
x02=x-x2
x03=x-x3
y01=y-y1
y02=y-y2
y03=y-y3

tmpcoef=(x01*y02+x02*y03+x03*y01-&
        &x03*y02-x02*y01-x01*y03)

coef(1)=(x02*y03-x03*y02)/tmpcoef
coef(2)=(x03*y01-x01*y03)/tmpcoef
coef(3)=(x01*y02-x02*y01)/tmpcoef

return
end subroutine tricoef
!-------------------------------------------------------------------------


! Calculate the coeficents of the interpolation of two nodes.
!-------------------------------------------------------------------------
subroutine bicoef(coef,x,y,x1,y1,x2,y2)

implicit none
real(kind=8) :: coef(4)
real(kind=8) :: x1,y1,x2,y2
real(kind=8) :: x,y,len1,len2


len1=sqrt((x-x1)**2+(y-y1)**2)
len2=sqrt((x2-x)**2+(y2-y1)**2)
coef(1)=len2/(len1+len2)
coef(2)=len1/(len1+len2)


return
end subroutine bicoef
!-------------------------------------------------------------------------




! This subroutine finds Frechet Derivative for inverse with regular 
! velocity grid.
!-----------------------------------------------------------------------
subroutine frech_regular(vel,path,ipath,fd,dx,dy,minx,maxx,miny,maxy,&
            &xnum,ynum)

implicit none
type(srstrct) :: path(maxpathnode)
type(tstrct), target :: vel(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: fd(maxvel)
real(kind=8) :: dx,dy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: x,y
real(kind=8) :: length
real(kind=8) :: a
integer :: ipath
integer :: xnum,ynum
integer :: i

allocate(prv(1:4))

fd=0d0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    length=(sqrt((path(i)%x-path(i+1)%x)**2&
          &+(path(i)%y-path(i+1)%y)**2))*radii&
          &*pi/1.8d2

    call localcood(x,y,vel,prv,dx,dy,minx,maxx,miny,maxy,xnum,ynum)

    a=(prv(4)%p%x-prv(1)%p%x)*(prv(4)%p%y-prv(1)%p%y)
    
    fd(prv(1)%p%num)=length*(prv(4)%p%x-x)*(prv(4)%p%y-y)/a&
                    &+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*(x-prv(1)%p%x)*(prv(4)%p%y-y)/a&
                    &+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*(prv(4)%p%x-x)*(y-prv(1)%p%y)/a&
                    &+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=length*(x-prv(1)%p%x)*(y-prv(1)%p%y)/a&
                    &+fd(prv(4)%p%num)

end do

deallocate(prv)

return
end subroutine frech_regular
!-----------------------------------------------------------------------




! This subroutine finds Frechet Derivative for inverse with 3d regular 
! velocity grid, in spherical coordinates.
!-----------------------------------------------------------------------
subroutine frech_regsph3d(vel,path,ipath,fd,dx,dy,dz,minx,maxx,miny,maxy,&
            &minz,maxz,xnum,ynum,znum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d), target :: vel(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: dx,dy,dz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: x,y,z
real(kind=8) :: length
real(kind=8) :: a
integer :: ipath
integer :: xnum,ynum,znum
integer :: i

allocate(prv(1:8))

fd=0d0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0
    
    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)

    call localcood3d(x,y,z,vel,prv,dx,dy,dz,minx,maxx,miny,maxy,&
                    &minz,maxz,xnum,ynum,znum)

    a=(prv(8)%p%x-prv(1)%p%x)*(prv(8)%p%y-prv(1)%p%y)*(prv(8)%p%z-prv(1)%p%z)
    
    fd(prv(1)%p%num)=length*(prv(8)%p%x-x)*(prv(8)%p%y-y)*(prv(8)%p%z-z)/a&
                    &+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*(x-prv(1)%p%x)*(prv(8)%p%y-y)*(prv(8)%p%z-z)/a&
                    &+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*(prv(8)%p%x-x)*(y-prv(1)%p%y)*(prv(8)%p%z-z)/a&
                    &+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=length*(x-prv(1)%p%x)*(y-prv(1)%p%y)*(prv(8)%p%z-z)/a&
                    &+fd(prv(4)%p%num)
    fd(prv(5)%p%num)=length*(prv(8)%p%x-x)*(prv(8)%p%y-y)*(z-prv(1)%p%z)/a&
                    &+fd(prv(5)%p%num)
    fd(prv(6)%p%num)=length*(x-prv(1)%p%x)*(prv(8)%p%y-y)*(z-prv(1)%p%z)/a&
                    &+fd(prv(6)%p%num)
    fd(prv(7)%p%num)=length*(prv(8)%p%x-x)*(y-prv(1)%p%y)*(z-prv(1)%p%z)/a&
                    &+fd(prv(7)%p%num)
    fd(prv(8)%p%num)=length*(x-prv(1)%p%x)*(y-prv(1)%p%y)*(z-prv(1)%p%z)/a&
                    &+fd(prv(8)%p%num)


end do

deallocate(prv)

return
end subroutine frech_regsph3d
!-----------------------------------------------------------------------




! This subroutine finds Frechet Derivative for inverse with 3d regular 
! velocity grid, in spherical coordinates.
!-----------------------------------------------------------------------
subroutine frech_regsph3d2(vel,path,ipath,fd,dx,dy,dz,minx,maxx,miny,maxy,&
            &minz,maxz,xnum,ynum,znum,topoxy,topoz,tpminx,tpmaxx,tpminy,&
            &tpmaxy,tpxnum,tpynum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d), target :: vel(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: dx,dy,dz
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: minz,maxz
real(kind=8) :: x,y,z
real(kind=8) :: length
real(kind=8) :: coflower(4),cofupper(4),zcof(2)
real(kind=8) :: tpminx,tpmaxx,tpminy,tpmaxy
real(kind=8) :: xtmp(4),ytmp(4),ztmp
integer :: lowerp(4),nlowerp
integer :: upperp(4),nupperp
integer :: ipath
integer :: xnum,ynum,znum
integer :: tpxnum,tpynum
integer :: i,j,k,updown

allocate(prv(1:8))

fd=0d0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0
    
    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)

    call localcood3d(x,y,z,vel,prv,dx,dy,dz,minx,maxx,miny,maxy,&
                    &minz,maxz,xnum,ynum,znum)

    
    ! Calculate the coefficiets of the nodes on the lower layer.
    ! Nodes in the air will not contribute.
    k=0
    do j=1,4
        call psurf(prv(j)%p%x,prv(j)%p%y,prv(j)%p%z,&
                  &topoxy,topoz,updown,tpminx,tpmaxx,&
                  &tpminy,tpmaxy,tpxnum,tpynum)
        if(updown .eq. 1)then
            k=k+1
            lowerp(k)=j
        end if
    end do
    nlowerp=k

    ! Covert the earth coordinate to earth surface local Cartesian coordinate.
    ! Because *coef subroutines use Cartesian coordinate.
    if(nlowerp .gt. 1)then
        xtmp(1)=0d0
        ytmp(1)=0d0
        do j=2,nlowerp
            call ear2loc(xtmp(j),ytmp(j),ztmp,prv(lowerp(j))%p%x,&
                       &prv(lowerp(j))%p%y,prv(lowerp(j))%p%z,&
                       &prv(lowerp(1))%p%x,prv(lowerp(1))%p%y,&
                       &prv(lowerp(1))%p%z)
        end do
    end if

    coflower=0d0
    select case(nlowerp)
    case(4)
        call quadcoef(coflower,x,y,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
    case(3)
        call tricoef(coflower,x,y,xtmp(1),ytmp(1),xtmp(2),&
                    &ytmp(2),xtmp(3),ytmp(3))
    case(2)
        call bicoef(coflower,x,y,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
    case(1)
        coflower(1)=1.0d0
    end select

    ! Calculate the coefficiets of the nodes on the upper layer.
    ! Nodes in the air will not contribute.
    k=0
    do j=5,8
        call psurf(prv(j)%p%x,prv(j)%p%y,prv(j)%p%z,&
                  &topoxy,topoz,updown,tpminx,tpmaxx,&
                  &tpminy,tpmaxy,tpxnum,tpynum)
        if(updown .eq. 1)then
            k=k+1
            upperp(k)=j
        end if
    end do
    nupperp=k


    ! Covert the earth coordinate to earth surface local Cartesian coordinate.
    ! Because *coef subroutines use Cartesian coordinate.
    ! xtmp, ytmp, ztmp are used the second time here.
    if(nupperp .gt. 1)then
        xtmp(1)=0d0
        ytmp(1)=0d0
        do j=2,nupperp
            call ear2loc(xtmp(j),ytmp(j),ztmp,prv(upperp(j))%p%x,&
                       &prv(upperp(j))%p%y,prv(upperp(j))%p%z,&
                       &prv(upperp(1))%p%x,prv(upperp(1))%p%y,&
                       &prv(upperp(1))%p%z)
        end do
    end if


    cofupper=0d0
    select case(nupperp)
    case(4)
        call quadcoef(cofupper,x,y,xtmp(1),ytmp(1),xtmp(4),ytmp(4))
    case(3)
        call tricoef(cofupper,x,y,xtmp(1),ytmp(1),xtmp(2),&
                    &ytmp(2),xtmp(3),ytmp(3))
    case(2)
        call bicoef(cofupper,x,y,xtmp(1),ytmp(1),xtmp(2),ytmp(2))
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
        zcof(1)=(prv(8)%p%z-z)/(prv(8)%p%z-prv(1)%p%z)
        zcof(2)=(z-prv(1)%p%z)/(prv(8)%p%z-prv(1)%p%z)
    end if
    
    fd(prv(1)%p%num)=length*zcof(1)*coflower(1)+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*zcof(1)*coflower(2)+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*zcof(1)*coflower(3)+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=length*zcof(1)*coflower(4)+fd(prv(4)%p%num)

    fd(prv(5)%p%num)=length*zcof(2)*cofupper(1)+fd(prv(5)%p%num)
    fd(prv(6)%p%num)=length*zcof(2)*cofupper(2)+fd(prv(6)%p%num)
    fd(prv(7)%p%num)=length*zcof(2)*cofupper(3)+fd(prv(7)%p%num)
    fd(prv(8)%p%num)=length*zcof(2)*cofupper(4)+fd(prv(8)%p%num)

end do

deallocate(prv)

return
end subroutine frech_regsph3d2
!-----------------------------------------------------------------------


! This subroutine finds Frechet Derivative for inverse with 3d regular 
! velocity grid, in spherical coordinates. Grid is non-uniform.
!-----------------------------------------------------------------------
subroutine frech_reg2sph3d(vel,path,ipath,fd,layer,nl)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d)  :: vel(maxgrid3d)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: layer(maxgrd1d,3)
real(kind=8) :: x,y,z
real(kind=8) :: length
real(kind=8) :: a
integer :: nl(3)
integer :: n(8)
integer :: ipath
integer :: i


fd=0d0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0
    
    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)

    call locatcood3d2(n,nl(1),layer(1:maxgrd1d,1),nl(2),layer(1:maxgrd1d,2)&
                    &,nl(3),layer(1:maxgrd1d,3),x,y,z)

    a=(vel(n(8))%x-vel(n(1))%x)*&
     &(vel(n(8))%y-vel(n(1))%y)*&
     &(vel(n(8))%z-vel(n(1))%z)

    
    fd(n(1))=length*(vel(n(8))%x-x)*(vel(n(8))%y-y)*(vel(n(8))%z-z)/a&
            &+fd(n(1))
    fd(n(2))=length*(x-vel(n(1))%x)*(vel(n(8))%y-y)*(vel(n(8))%z-z)/a&
            &+fd(n(2))
    fd(n(3))=length*(vel(n(8))%x-x)*(y-vel(n(1))%y)*(vel(n(8))%z-z)/a&
            &+fd(n(3))
    fd(n(4))=length*(x-vel(n(1))%x)*(y-vel(n(1))%y)*(vel(n(8))%z-z)/a&
            &+fd(n(4))
    fd(n(5))=length*(vel(n(8))%x-x)*(vel(n(8))%y-y)*(z-vel(n(1))%z)/a&
            &+fd(n(5))
    fd(n(6))=length*(x-vel(n(1))%x)*(vel(n(8))%y-y)*(z-vel(n(1))%z)/a&
            &+fd(n(6))
    fd(n(7))=length*(vel(n(8))%x-x)*(y-vel(n(1))%y)*(z-vel(n(1))%z)/a&
            &+fd(n(7))
    fd(n(8))=length*(x-vel(n(1))%x)*(y-vel(n(1))%y)*(z-vel(n(1))%z)/a&
            &+fd(n(8))


end do


return
end subroutine frech_reg2sph3d
!-----------------------------------------------------------------------


! This subroutine finds Frechet Derivative for inverse with 3d regular 
! velocity grid, in spherical coordinates.
!-----------------------------------------------------------------------
subroutine frech_reg2sph3d2(vel,path,ipath,fd,layer,nl,&
            &topoxy,topoz,tpminx,tpmaxx,tpminy,tpmaxy,tpxnum,tpynum)

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
integer :: lowerp(4),nlowerp
integer :: upperp(4),nupperp
integer :: nl(3)
integer :: n(8)
integer :: ipath
integer :: tpxnum,tpynum
integer :: i,j,k,updown


fd=0d0

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

    
    fd(n(1))=length*zcof(1)*coflower(1)+fd(n(1))
    fd(n(2))=length*zcof(1)*coflower(2)+fd(n(2))
    fd(n(3))=length*zcof(1)*coflower(3)+fd(n(3))
    fd(n(4))=length*zcof(1)*coflower(4)+fd(n(4))

    fd(n(5))=length*zcof(2)*cofupper(1)+fd(n(5))
    fd(n(6))=length*zcof(2)*cofupper(2)+fd(n(6))
    fd(n(7))=length*zcof(2)*cofupper(3)+fd(n(7))
    fd(n(8))=length*zcof(2)*cofupper(4)+fd(n(8))

end do


return
end subroutine frech_reg2sph3d2
!-----------------------------------------------------------------------



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
fd(1+(evnid-1)*4)=1.0d0
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

    
    fd(n(1)+ttlevn*4)=length*zcof(1)*coflower(1)+fd(n(1)+ttlevn*4)
    fd(n(2)+ttlevn*4)=length*zcof(1)*coflower(2)+fd(n(2)+ttlevn*4)
    fd(n(3)+ttlevn*4)=length*zcof(1)*coflower(3)+fd(n(3)+ttlevn*4)
    fd(n(4)+ttlevn*4)=length*zcof(1)*coflower(4)+fd(n(4)+ttlevn*4)

    fd(n(5)+ttlevn*4)=length*zcof(2)*cofupper(1)+fd(n(5)+ttlevn*4)
    fd(n(6)+ttlevn*4)=length*zcof(2)*cofupper(2)+fd(n(6)+ttlevn*4)
    fd(n(7)+ttlevn*4)=length*zcof(2)*cofupper(3)+fd(n(7)+ttlevn*4)
    fd(n(8)+ttlevn*4)=length*zcof(2)*cofupper(4)+fd(n(8)+ttlevn*4)

end do


return
end subroutine frech_reg2sph3d2rl
!-----------------------------------------------------------------------



! This subroutine finds Frechet Derivative for inverse with Delaunay
! triangle velocity grid.
!-----------------------------------------------------------------------
subroutine frech_tri(vel,path,ipath,fd,tri,trinum)

implicit none
type(srstrct) :: path(maxpathnode)
type(tstrct), target :: vel(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: fd(maxvel)
real(kind=8) :: x,y
real(kind=8) :: length
real(kind=8) :: area
real(kind=8) :: area1,area2,area3
integer :: tri(maxtri,3)
integer :: trinum
integer :: ipath
integer :: i

allocate(prv(1:3))

fd=0d0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    length=(sqrt((path(i)%x-path(i+1)%x)**2&
          &+(path(i)%y-path(i+1)%y)**2))*radii&
          &*pi/1.8d2

    call baryc_locat(x,y,vel,tri,trinum,prv)
    area=abs(prv(1)%p%x*prv(2)%p%y - prv(2)%p%x*prv(1)%p%y&
            &+ prv(2)%p%x*prv(3)%p%y - prv(3)%p%x*prv(2)%p%y&
            &+ prv(3)%p%x*prv(1)%p%y - prv(1)%p%x*prv(3)%p%y)
    
    area1=abs(x*prv(2)%p%y - prv(2)%p%x*y&
            &+ prv(2)%p%x*prv(3)%p%y - prv(3)%p%x*prv(2)%p%y&
            &+ prv(3)%p%x*y - x*prv(3)%p%y)

    area2=abs(prv(1)%p%x*y - x*prv(1)%p%y&
            &+ x*prv(3)%p%y - prv(3)%p%x*y&
            &+ prv(3)%p%x*prv(1)%p%y - prv(1)%p%x*prv(3)%p%y)
    
    area3=area-area1-area2
    fd(prv(1)%p%num)=length*area1/area+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*area2/area+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*area3/area+fd(prv(3)%p%num)

end do

deallocate(prv)

return
end subroutine frech_tri
!-----------------------------------------------------------------




! This subroutine finds Frechet Derivative for inverse with 3D
! tetrahedro velocity grid, in sphirical coordinates.
!-----------------------------------------------------------------
subroutine frech_tethsph(vel,path,ipath,fd,teth,tethnbr,tethnum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d), target :: vel(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: x,y,z
real(kind=8) :: length
real(kind=8) :: vol
real(kind=8) :: vol1,vol2,vol3,vol4
integer :: teth(maxteth,4)
integer :: tethnbr(maxteth,4)
integer :: tethnum
integer :: ipath
integer :: i
integer :: itet

allocate(prv(1:4))

fd=0d0

do i=1,ipath-1

    ! Slowness is chosen at the center of path segment.
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0

    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)

    call dllocat3d(x,y,z,vel,teth,tethnbr,tethnum,prv,itet)


    call teth_volsph(vol,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                        &prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                        &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z,&
                        &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol1,x,y,z,prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                               &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol2,x,y,z,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                               &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol3,x,y,z,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                               &prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    vol4=vol-vol1-vol2-vol3

    
    fd(prv(1)%p%num)=length*vol1/vol+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*vol2/vol+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*vol3/vol+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=length*vol4/vol+fd(prv(4)%p%num)

end do

deallocate(prv)

return
end subroutine frech_tethsph
!-----------------------------------------------------------------



! This subroutine finds Frechet Derivative for inverse with 3D
! tetrahedro velocity grid.
!-----------------------------------------------------------------
subroutine frech_tethsph2(vel,path,ipath,fd,teth,tethnbr,tethnum,&
                      &topoxy,topoz,tpminx,tpmaxx,&
                      &tpminy,tpmaxy,tpxnum,tpynum)

implicit none
type(srstrct3d) :: path(maxpathnode3d)
type(tstrct3d), target :: vel(maxgrid3d)
type(pstrct3d), pointer :: prv(:)
type(tstrct) :: topoxy(maxgrid)
real(kind=8) :: topoz(maxgrid)
real(kind=8) :: fd(maxvel3d)
real(kind=8) :: x,y,z
real(kind=8) :: length
real(kind=8) :: vol
real(kind=8) :: vol1,vol2,vol3,vol4
real(kind=8) :: tpminx,tpmaxx,tpminy,tpmaxy
integer :: vtpos(4)
integer :: teth(maxteth,4)
integer :: tethnbr(maxteth,4)
integer :: tethnum
integer :: tpxnum,tpynum
integer :: ipath
integer :: i,j
integer :: itet
integer :: updown

allocate(prv(1:4))

fd=0d0

do i=1,ipath-1

    ! Slowness is chosen at the center of path segment.
    x=(path(i)%x+path(i+1)%x)/2.0d0
    y=(path(i)%y+path(i+1)%y)/2.0d0
    z=(path(i)%z+path(i+1)%z)/2.0d0

    call sphdist(path(i)%x,path(i)%y,path(i)%z,&
                &path(i+1)%x,path(i+1)%y,path(i+1)%z,&
                &length)
                
    call dllocat3d(x,y,z,vel,teth,tethnbr,tethnum,prv,itet)

    do j=1,4
        call psurf(prv(j)%p%x,prv(j)%p%y,prv(j)%p%z,&
                 &topoxy,topoz,updown,tpminx,tpmaxx,&
                 &tpminy,tpmaxy,tpxnum,tpynum)
        if((updown .eq. 1) .or. (updown .eq. 0))then
            vtpos(j)=1
        else if(updown .eq. -1)then
            vtpos(j)=0
        else
            write(*,*)"Err in strct.f90: frech_teth2: psurf! &
                      &parameter updown = ",updown
        end if
    end do

    call teth_volsph(vol1,x,y,z,prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                               &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol2,x,y,z,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                               &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol3,x,y,z,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                               &prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                               &prv(4)%p%x,prv(4)%p%y,prv(4)%p%z)

    call teth_volsph(vol4,x,y,z,prv(1)%p%x,prv(1)%p%y,prv(1)%p%z,&
                               &prv(2)%p%x,prv(2)%p%y,prv(2)%p%z,&
                               &prv(3)%p%x,prv(3)%p%y,prv(3)%p%z)


    vol=vol1*dble(vtpos(1))+vol2*dble(vtpos(2))&
      &+vol3*dble(vtpos(3))+vol4*dble(vtpos(4))

    
    fd(prv(1)%p%num)=dble(vtpos(1))*length*vol1/vol+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=dble(vtpos(2))*length*vol2/vol+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=dble(vtpos(3))*length*vol3/vol+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=dble(vtpos(4))*length*vol4/vol+fd(prv(4)%p%num)

end do

deallocate(prv)

return
end subroutine frech_tethsph2
!-----------------------------------------------------------------



! This subroutine determines 4 points in on a plane or on a
! tetrahedro.
!-----------------------------------------------------------------
subroutine teth_vol(vol,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

implicit none
real(kind=8) :: x1,y1,z1
real(kind=8) :: x2,y2,z2
real(kind=8) :: x3,y3,z3
real(kind=8) :: x4,y4,z4
real(kind=8) :: vol

vol=abs((x1-x4)*(y2-y4)*(z3-z4) + (y1-y4)*(z2-z4)*(x3-x4)&
      &+(z1-z4)*(x2-x4)*(y3-y4) - (x3-x4)*(y2-y4)*(z1-z4)&
      &-(y3-y4)*(z2-z4)*(x1-x4) - (z3-z4)*(x2-x4)*(y1-y4))/6.0d0


return
end subroutine teth_vol
!-----------------------------------------------------------------




! This subroutine determines if 4 points on a plane or not
!-----------------------------------------------------------------
subroutine teth_volsph(vol,lon1,lan1,r1,lon2,lan2,r2,&
                      &lon3,lan3,r3,lon4,lan4,r4)

implicit none
real(kind=8) :: lon1,lon2,lon3,lon4
real(kind=8) :: lan1,lan2,lan3,lan4
real(kind=8) :: r1,r2,r3,r4
real(kind=8) :: x1,x2,x3,x4
real(kind=8) :: y1,y2,y3,y4
real(kind=8) :: z1,z2,z3,z4
real(kind=8) :: vol

call sph2xyz(x1,y1,z1,lon1,lan1,r1)
call sph2xyz(x2,y2,z2,lon2,lan2,r2)
call sph2xyz(x3,y3,z3,lon3,lan3,r3)
call sph2xyz(x4,y4,z4,lon4,lan4,r4)

vol=abs((x1-x4)*(y2-y4)*(z3-z4)&
      &+(y1-y4)*(z2-z4)*(x3-x4)&
      &+(z1-z4)*(x2-x4)*(y3-y4)&
      &-(x3-x4)*(y2-y4)*(z1-z4)&
      &-(y3-y4)*(z2-z4)*(x1-x4)&
      &-(z3-z4)*(x2-x4)*(y1-y4))/6.0d0


return
end subroutine teth_volsph
!-----------------------------------------------------------------


! This subroutine add nodes in Delaunay triangles with high ray weights.
! (With triangles which have at least two high ray weights vertexes.)
! Input: nodenum (nodes's numbers which ray weights are high)
!        num (the total num of 'nodenum')
!        weight (ray weight of velocity nodes)
!        tri (Delaunay triangles), 
!        trinum (numbers of D tri), 
!        velnode (velocity nodes),
! Output: x (x coordinates of new nodes which will be added)
!         y (y coordinates of new nodes which will be added)
!         i (the total number of new nodes which will be added)
!-----------------------------------------------------------------------
subroutine addnodes(x,y,i,nodenum,num,weight,velnode,tri,trinum)

implicit none
type(tstrct), target :: velnode(maxgrid)
type(gnstrct) :: weight(maxvel)
real(kind=8) :: x(maxtri),y(maxtri)
real(kind=8) :: weightsum
integer :: tri(maxtri,3)
integer :: nodenum(maxvel)
integer :: it(3)
integer :: itri,trinum
integer :: ind,num
integer :: i,j


i=0
do itri=1,trinum
    j=0
    it=0
    do ind=1,num
        if(nodenum(ind) .eq. (tri(itri,1)+1))then
            j=j+1
            it(1)=it(1)+1
            if(it(1) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (tri(itri,2)+1))then
            j=j+1
            it(2)=it(2)+1
            if(it(2) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (tri(itri,3)+1))then
            j=j+1
            it(3)=it(3)+1
            if(it(3) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        end if
    end do
    if(j .ge. 2)then
        i=i+1
        weightsum=weight(tri(itri,1)+1)%val+&
                 &weight(tri(itri,2)+1)%val+&
                 &weight(tri(itri,3)+1)%val
        x(i)=(weight(tri(itri,1)+1)%val*velnode(tri(itri,1)+1)%x&
            &+weight(tri(itri,2)+1)%val*velnode(tri(itri,2)+1)%x&
            &+weight(tri(itri,3)+1)%val*velnode(tri(itri,3)+1)%x&
            &)/weightsum
        y(i)=(weight(tri(itri,1)+1)%val*velnode(tri(itri,1)+1)%y&
            &+weight(tri(itri,2)+1)%val*velnode(tri(itri,2)+1)%y&
            &+weight(tri(itri,3)+1)%val*velnode(tri(itri,3)+1)%y&
            &)/weightsum
    end if
end do


return
end subroutine addnodes
!-----------------------------------------------------------------



! This subroutine add nodes in tetrahedros with high ray weights.
! (With tetrahedros which have at least two high ray weights vertexes.)
! Input: nodenum (numbers of nodes which ray weights are high)
!        num (the total num of 'nodenum')
!        weight (ray weight of velocity nodes)
!        teth (Tetrahedro), 
!        trinum (numbers of Teths), 
!        velnode (adaptive velocity nodes),
!        dthrshd (minimum distance between two new nodes)
! Output: x (x coordinates of new nodes which will be added)
!         y (y coordinates of new nodes which will be added)
!         z (z coordinates of new nodes which will be added)
!         i (the total number of new nodes which will be added)
!-----------------------------------------------------------------------
subroutine addnodes3d(x,y,z,i,nodenum,num,weight,velnode,teth,tethnum,&
                     &dthrshd)

implicit none
type(tstrct3d), target :: velnode(maxgrid3d)
type(gnstrct) :: weight(maxvel3d)
real(kind=8) :: x(maxteth),y(maxteth),z(maxteth)
real(kind=8) :: weightsum
real(kind=8) :: dthrshd,dist
real(kind=8) :: xtemp,ytemp,ztemp
integer :: teth(maxteth,4)
integer :: nodenum(maxvel3d)
integer :: it(4)
integer :: iteth,tethnum
integer :: ind,num
integer :: i,j
integer :: status1,idst

i=0
do iteth=1,tethnum
    j=0
    it=0
    do ind=1,num
        if(nodenum(ind) .eq. (teth(iteth,1)+1))then
            j=j+1
            it(1)=it(1)+1
            if(it(1) .ge. 2)then
                write(*,*)"Error in addnode3d, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (teth(iteth,2)+1))then
            j=j+1
            it(2)=it(2)+1
            if(it(2) .ge. 2)then
                write(*,*)"Error in addnode3d, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (teth(iteth,3)+1))then
            j=j+1
            it(3)=it(3)+1
            if(it(3) .ge. 2)then
                write(*,*)"Error in addnode3d, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (teth(iteth,4)+1))then
            j=j+1
            it(4)=it(4)+1
            if(it(4) .ge. 2)then
                write(*,*)"Error in addnode3d, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        end if
    end do
    if(j .ge. 2)then
        weightsum=weight(teth(iteth,1)+1)%val+&
                 &weight(teth(iteth,2)+1)%val+&
                 &weight(teth(iteth,3)+1)%val+&
                 &weight(teth(iteth,4)+1)%val
        xtemp=(weight(teth(iteth,1)+1)%val*velnode(teth(iteth,1)+1)%x&
            &+weight(teth(iteth,2)+1)%val*velnode(teth(iteth,2)+1)%x&
            &+weight(teth(iteth,3)+1)%val*velnode(teth(iteth,3)+1)%x&
            &+weight(teth(iteth,4)+1)%val*velnode(teth(iteth,4)+1)%x&
            &)/weightsum
        ytemp=(weight(teth(iteth,1)+1)%val*velnode(teth(iteth,1)+1)%y&
            &+weight(teth(iteth,2)+1)%val*velnode(teth(iteth,2)+1)%y&
            &+weight(teth(iteth,3)+1)%val*velnode(teth(iteth,3)+1)%y&
            &+weight(teth(iteth,4)+1)%val*velnode(teth(iteth,4)+1)%y&
            &)/weightsum
        ztemp=(weight(teth(iteth,1)+1)%val*velnode(teth(iteth,1)+1)%z&
            &+weight(teth(iteth,2)+1)%val*velnode(teth(iteth,2)+1)%z&
            &+weight(teth(iteth,3)+1)%val*velnode(teth(iteth,3)+1)%z&
            &+weight(teth(iteth,4)+1)%val*velnode(teth(iteth,4)+1)%z&
            &)/weightsum

        ! The new nodes whose distances with previous nodes are greater
        ! than dthrshd will be accepted.
        if(i .ge. 1)then
            status1=1
            do idst=1,i
                call sphdist(xtemp,ytemp,ztemp,x(idst),y(idst),z(idst),dist)
                if(dist .lt. dthrshd)then
                    status1=0
                    exit
                end if
            end do

            if(status1 .eq. 1)then
                i=i+1
                x(i)=xtemp
                y(i)=ytemp
                z(i)=ztemp
            elseif(status1 .eq. 0)then
                write(*,*)"Distance between two node is too small, node1, node2, dist",&
                        &xtemp,ytemp,ztemp,x(idst),y(idst),z(idst),dist
            else
                write(*,*)"Error in subroutine addnodes3 of strct.f90,&
                         &status1=",status1
            end if
        elseif(i .eq. 0)then
            i=i+1
            x(i)=xtemp
            y(i)=ytemp
            z(i)=ztemp
        else
            write(*,*)"Error in subroutine addnodes3 of strct.f90,&
                     &i=",i
        end if
    end if
end do


return
end subroutine addnodes3d
!-----------------------------------------------------------------


! Calculate the distance between two points in spherical coordinates.
!------------------------------------------------------------------
subroutine sphdist(lon1,lan1,z1,lon2,lan2,z2,dist)

implicit none
real(kind=8) :: lan1,lan2
real(kind=8) :: lon1,lon2
real(kind=8) :: z1,z2
real(kind=8) :: tmp1,tmp2
real(kind=8) :: theta
real(kind=8) :: r1,r2
real(kind=8) :: dist


tmp1=(sin((lan1-lan2)*pi/3.6d2))**2
tmp2=((sin((lon1-lon2)*pi/3.6d2))**2)*cos(lan1*pi/1.8d2)*cos(lan2*pi/1.8d2)

theta=2.0d0*asin(sqrt(tmp1+tmp2))
r1=radii-z1
r2=radii-z2

dist=sqrt(r1*r1+r2*r2-2.0d0*r1*r2*cos(theta))

return
end subroutine sphdist

!-------------------------------------------------------------------------



! Covert sphirical coordinate to Cartesian coordinate.
!------------------------------------------------------------------
subroutine sph2xyz(x,y,z,theta,phi,r)

implicit none
real(kind=8) :: theta,phi,r
real(kind=8) :: x,y,z


x=r*cos(deg2rad(theta))*cos(deg2rad(phi))
y=r*cos(deg2rad(theta))*sin(deg2rad(phi))
z=r*sin(deg2rad(theta))

return
end subroutine sph2xyz

!-------------------------------------------------------------------------


! Covert Cartesian coordinate to sphirical coordinate.
!------------------------------------------------------------------
subroutine xyz2sph(theta,phi,r,x,y,z)

implicit none
real(kind=8) :: theta,phi,r
real(kind=8) :: x,y,z

r=sqrt(x**2+y**2+z**2)
theta=rad2deg(asin(z/r))
if(x .ge. 0d0 .and. y .ge. 0d0)then     
    phi=rad2deg(atan(y/x))          
elseif(x .lt. 0d0 .and. y .ge. 0d0)then 
    phi=1.8d0-rad2deg(atan(-y/x))   
elseif(x .lt. 0d0 .and. y .lt. 0d0)then 
    phi=1.8d0+rad2deg(atan(y/x))    
else                                
    phi=3.6d0+rad2deg(atan(y/x))    
end if

return
end subroutine xyz2sph

!-------------------------------------------------------------------------


! Covert earth coordinate to Cartesian coordinate.
!------------------------------------------------------------------
subroutine ear2xyz(x,y,z,lon,lan,dep)

implicit none
real(kind=8) :: lon,lan,dep
real(kind=8) :: x,y,z,r

r=radii-dep
x=r*sin(deg2rad(lan))*cos(deg2rad(lon))
y=r*sin(deg2rad(lan))*sin(deg2rad(lon))
z=r*cos(deg2rad(lan))

return

end subroutine ear2xyz
!-----------------------------------------------------------------


! Covert Cartesian coordinate to earth coordinate.
!------------------------------------------------------------------
subroutine xyz2ear(lon,lan,dep,x,y,z)

implicit none
real(kind=8) :: lon,lan,dep
real(kind=8) :: x,y,z,r

r=sqrt(x**2+y**2+z**2)
dep=radii-r
lan=9.0d0-rad2deg(acos(z/r))

if(x .ge. 0d0)then
    if(y .ge. 0d0)then
        lon=rad2deg(atan(y/x))
    else
        lon=1.8d2-rad2deg(atan(y/x))
    end if
else
    if(y .ge. 0d0)then
        lon=1.8d2-rad2deg(atan(y/x))
    else
        lon=3.6d2+rad2deg(atan(y/x))
    end if
end if

return
end subroutine xyz2ear

!-------------------------------------------------------------------------


! Covert earth coordinate to earth surface local Cartesian coordinate.
!------------------------------------------------------------------
subroutine ear2loc(x,y,z,lon,lan,dep,lon0,lan0,dep0)

implicit none
real(kind=8) :: lon,lan,dep,lon0,lan0,dep0
real(kind=8) :: x,y,z,r

r=radii-dep
x=r*cos(deg2rad(lan))*deg2rad(lon-lon0)
y=r*deg2rad(lan-lan0)
z=dep-dep0

return

end subroutine ear2loc
!-----------------------------------------------------------------


! Covert Cartesian coordinate to earth coordinate.
!------------------------------------------------------------------
subroutine loc2ear(lon,lan,dep,lon0,lan0,dep0,x,y,z)

implicit none
real(kind=8) :: lon,lan,dep,lon0,lan0,dep0
real(kind=8) :: x,y,z,r

dep=z+dep0
r=radii-dep
lan=rad2deg(y/r)+lan0
lon=rad2deg(x/(r*cos(deg2rad(lan))))+lon0

return
end subroutine loc2ear

!-------------------------------------------------------------------------


! Calculate the distance between two points in Cartesian coordinates.
!------------------------------------------------------------------
subroutine xyzdist(x1,y1,z1,x2,y2,z2,dist)

implicit none
real(kind=8) :: x1,x2
real(kind=8) :: y1,y2
real(kind=8) :: z1,z2
real(kind=8) :: dist


dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)


end subroutine xyzdist

!-------------------------------------------------------------------------


!Determine whehter a point(x,y,z) is above the surface or not.
!updown=1,up or on surface; -1,down
!--------------------------------------------------------------------------
subroutine psurf(x,y,z,surfxy,surfz,updown,sminx,smaxx,sminy,smaxy,sxnum,synum)


implicit none
type(tstrct), target :: surfxy(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: surfz(maxgrid)
real(kind=8) :: x,y,z
real(kind=8) :: tx(3),ty(3),tz(3)
real(kind=8) :: sminx,smaxx
real(kind=8) :: sminy,smaxy
real(kind=8) :: sdx,sdy
real(kind=8) :: a
real(kind=8) :: normz
real(kind=8) :: cof
real(kind=8) :: fsurf
integer :: sxnum,synum
integer :: updown

allocate(prv(1:4))
sdx=(smaxx-sminx)/dble(sxnum-1)
sdy=(smaxy-sminy)/dble(synum-1)
call localcood(x,y,surfxy,prv,sdx,sdy,sminx,smaxx,sminy,smaxy,sxnum,synum)


a=(prv(4)%p%y-prv(1)%p%y)/(prv(4)%p%x-prv(1)%p%x)

tx(1)=prv(1)%p%x
tx(3)=prv(4)%p%x
ty(1)=prv(1)%p%y
ty(3)=prv(4)%p%y
tz(1)=surfz(prv(1)%p%num)
tz(3)=surfz(prv(4)%p%num)
if( (y-a*x) .ge. 0d0)then
    tx(2)=prv(3)%p%x
    ty(2)=prv(3)%p%y
    tz(2)=surfz(prv(3)%p%num)
else
    tx(2)=prv(2)%p%x
    ty(2)=prv(2)%p%y
    tz(2)=surfz(prv(2)%p%num)
end if


normz=(tx(2)-tx(1))*(ty(3)-ty(1))-(tx(3)-tx(1))*(ty(2)-ty(1))

if(normz .ge. 0d0)then
    cof=1d0
else
    cof=-1d0
end if

fsurf=cof*((x-tx(1))*(ty(2)-ty(1))*(tz(3)-tz(1))&
         &+(y-ty(1))*(tz(2)-tz(1))*(tx(3)-tx(1))&
         &+(z-tz(1))*(tx(2)-tx(1))*(ty(3)-ty(1))&
         &-(x-tx(1))*(ty(3)-ty(1))*(tz(2)-tz(1))&
         &-(y-ty(1))*(tz(3)-tz(1))*(tx(2)-tx(1))&
         &-(z-tz(1))*(tx(3)-tx(1))*(ty(2)-ty(1)))

if(fsurf .ge. 0d0)then
    updown=1
else
    updown=-1
end if

deallocate(prv)

return
end subroutine psurf
!------------------------------------------------------------------------



!Determine the elevation of a point(x,y,z) on the surface.
!--------------------------------------------------------------------------
subroutine surfele(x,y,z,surfxy,surfz,sminx,smaxx,sminy,smaxy,sxnum,synum)


implicit none
type(tstrct), target :: surfxy(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: surfz(maxgrid)
real(kind=8) :: x,y,z
real(kind=8) :: tx(3),ty(3),tz(3)
real(kind=8) :: sminx,smaxx
real(kind=8) :: sminy,smaxy
real(kind=8) :: sdx,sdy
real(kind=8) :: a
integer :: sxnum,synum


allocate(prv(1:4))
sdx=(smaxx-sminx)/dble(sxnum-1)
sdy=(smaxy-sminy)/dble(synum-1)
call localcood(x,y,surfxy,prv,sdx,sdy,sminx,smaxx,sminy,smaxy,sxnum,synum)


a=(prv(4)%p%y-prv(1)%p%y)/(prv(4)%p%x-prv(1)%p%x)

tx(1)=prv(1)%p%x
tx(3)=prv(4)%p%x
ty(1)=prv(1)%p%y
ty(3)=prv(4)%p%y
tz(1)=surfz(prv(1)%p%num)
tz(3)=surfz(prv(4)%p%num)
if( (y-a*x) .ge. 0d0)then
    tx(2)=prv(3)%p%x
    ty(2)=prv(3)%p%y
    tz(2)=surfz(prv(3)%p%num)
else
    tx(2)=prv(2)%p%x
    ty(2)=prv(2)%p%y
    tz(2)=surfz(prv(2)%p%num)
end if


z=trisurfintp(tz(1),tz(2),tz(3),tx(1),ty(1),tx(2),ty(2),tx(3),ty(3),x,y)

deallocate(prv)

return
end subroutine surfele
!------------------------------------------------------------------------



! Merge a list of files into a single file.
! pfile: input file list; nf: number of input files;
! sfile: output single file
! ds=0: delete original files; ds=1: keep files
!-------------------------------------------------------------------------
subroutine mergefile(pfile,nf,sfile,ds)

implicit none
integer :: i,nf
integer :: ds
integer :: status1
character(len=70) :: pfile(maxsource3d)
character(len=*) :: sfile
character(len=200) :: tmpchar

open(201,file=sfile,status='replace')
do i=1,nf
    open(301,file=trim(pfile(i)),status='old')
    do while(.true.)
        read(301,'(a)',iostat=status1)tmpchar
        if(status1/=0)exit
        write(201,*)trim(tmpchar)
    end do
    if(ds .eq. 0)then
        close(301,status='delete')
    else if(ds .eq. 1)then
        close(301)
    else
        write(*,*)"error in mergefile of strct.f90!"
        write(*,*)"ds = ",ds
        stop
    end if
end do
close(201)

return
end subroutine mergefile
!-------------------------------------------------------------------------



! Merge a list of binary files into a single file.
! pfile: input file list; nf: number of input files;
! sfile: output single file
! ds=0: delete original files; ds=1: keep files
!-------------------------------------------------------------------------
subroutine mergebf(pfile,nf,reclength,sfile,ds)

implicit none
integer :: i,j,k,nf
integer :: ds,reclength
integer :: status1
character(len=70) :: pfile(maxsource3d)
character(len=*) :: sfile
character(len=reclength) :: tmpchar

open(201,file=sfile,status='replace',form='unformatted',&
    &access='direct',recl=reclength)
j=0
do i=1,nf
    open(301,file=trim(pfile(i)),status='old',form='unformatted',&
        &access='direct',recl=reclength)
    k=0
    do while(.true.)
        k=k+1
        read(301,rec=k,iostat=status1)tmpchar
        if(status1/=0)exit
        j=j+1
        write(201,rec=j)tmpchar
    end do
    if(ds .eq. 0)then
        close(301,status='delete')
    else if(ds .eq. 1)then
        close(301)
    else
        write(*,*)"error in mergebf of strct.f90!"
        write(*,*)"ds = ",ds
        stop
    end if
end do
close(201)

return
end subroutine mergebf
!-------------------------------------------------------------------------



! Calculate the direction vector of a line.
! coord=1: Cartesian coordinate; coord=2: spherical coordinate.
! When use spherical coordinate: x is longitude, y is lantitude.
!------------------------------------------------------------------------
subroutine directionvector(x1,y1,z1,x2,y2,z2,xv,yv,zv,coord)

implicit none
real(kind=8) :: x1,y1,z1
real(kind=8) :: x2,y2,z2
real(kind=8) :: xv,yv,zv
real(kind=8) :: ds
integer :: coord

if(coord .eq. 1)then
    ds=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
else if(coord .eq. 2)then
    call sphdist(x1,y1,z1,x2,y2,z2,ds)
else
    write(*,*)"error in directionvector of strct.f90!"
    write(*,*)"coord = ",coord
end if

xv=(x2-x1)/ds
yv=(y2-y1)/ds
zv=(z2-z1)/ds

return
end subroutine directionvector
!------------------------------------------------------------------------


! Calculate velocity of an point from 1D velocity model.
!-------------------------------------------------------------------------
subroutine vfrom1d(v,depth,vinp,vnum)

implicit none
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: depth
real(kind=8) :: v
integer :: iv,vnum
integer :: status1

status1=0
do iv=1,vnum
    if((depth .ge. vinp(iv,1)) .and. (depth .lt. vinp(iv,2)))then
        v=vinp(iv,3)
        status1=status1+1
        exit
    end if
end do

if(status1/=1)then
    write(*,*)"Error in vfrom1d in strct.f90!"
    write(*,*)"Status1=",status1,"V=",v
    write(*,*)"Depth=",depth
    stop
end if


return
end subroutine vfrom1d

!-------------------------------------------------------------------------


! Calculate velocity of an point from 1D velocity model.
!-------------------------------------------------------------------------
subroutine vfrom1dct(v,depth,vinp,vnum)

implicit none
real(kind=8) :: vinp(maxgrd1d,3)
real(kind=8) :: depth
real(kind=8) :: v
integer :: iv,vnum
integer :: status1

status1=0
do iv=1,vnum
    if((depth .ge. vinp(iv,1)) .and. (depth .lt. vinp(iv,2)))then
        v=(vinp(iv,3)*(vinp(iv,2)-depth)+&
          &vinp(iv+1,3)*(depth-vinp(iv,1)))/&
          &(vinp(iv,2)-vinp(iv,1))
        status1=status1+1
        exit
    end if
end do

if(status1/=1)then
    write(*,*)"Error in vfrom1dct in strct.f90!"
    write(*,*)"Status1=",status1,"V=",v
    write(*,*)"Depth=",depth
    stop
end if


return
end subroutine vfrom1dct
!-------------------------------------------------------------------------


! Calculate determinant of 3x3 matrix
!-------------------------------------------------------------------------
subroutine det3(adet,a11,a12,a13,a21,a22,a23,a31,a32,a33)

implicit none
real(kind=8) :: a11,a12,a13,a21,a22,a23,a31,a32,a33
real(kind=8) :: adet

adet=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31)

return
end subroutine det3

!-------------------------------------------------------------------------------


! Calculate determinant of 4x4 matrix
!-------------------------------------------------------------------------
subroutine det4(a,adet)

implicit none
real(kind=8) :: a(4,4)
real(kind=8) :: m(4)
real(kind=8) :: adet


call det3(m(1),a(2,2),a(2,3),a(2,4),&
              &a(3,2),a(3,3),a(3,4),&
              &a(4,2),a(4,3),a(4,4))
call det3(m(2),a(1,2),a(1,3),a(1,4),&
              &a(3,2),a(3,3),a(3,4),&
              &a(4,2),a(4,3),a(4,4))
call det3(m(3),a(1,2),a(1,3),a(1,4),&
              &a(2,2),a(2,3),a(2,4),&
              &a(4,2),a(4,3),a(4,4))
call det3(m(4),a(1,2),a(1,3),a(1,4),&
              &a(2,2),a(2,3),a(2,4),&
              &a(3,2),a(3,3),a(3,4))

adet=a(1,1)*m(1)-a(2,1)*m(2)+a(3,1)*m(3)-a(4,1)*m(4)

return
end subroutine det4

!-------------------------------------------------------------------------------

! This subroutine locates the point (x,y,z) in tethahedrons.
!----------------------------------------------------------------------
subroutine dllocat3d(x,y,z,node,teth,tethnbr,tethnum,pteth,itet)

implicit none
type(tstrct3d), target :: node(maxgrid3d)
type(pstrct3d), pointer :: pteth(:)
real(kind=8) :: x,y,z
real(kind=8) :: x1,x2,x3,xp
real(kind=8) :: y1,y2,y3,yp
real(kind=8) :: z1,z2,z3,zp
real(kind=8) :: vol,vol1,vol2,vol3,vol4
integer :: teth(maxteth,4)
integer :: tethnbr(maxteth,4)
integer :: tethnum
integer :: it,i,iloop
integer :: itet
integer :: j

allocate(pteth(1:4))
do i=1,4
    nullify(pteth(i)%p)
end do


j=0
it=1
iloop=0
do while((iloop .eq. 0) .and. (j .le. tethnum))
    j=j+1 
    !write(*,*)"line1956 in strct.f90:", it
    x1=(node(teth(it,1)+1)%x-node(teth(it,4)+1)%x)*1.0d2
    y1=(node(teth(it,1)+1)%y-node(teth(it,4)+1)%y)*1.1d2
    z1=node(teth(it,1)+1)%z-node(teth(it,4)+1)%z
    x2=(node(teth(it,2)+1)%x-node(teth(it,4)+1)%x)*1.0d2
    y2=(node(teth(it,2)+1)%y-node(teth(it,4)+1)%y)*1.1d2
    z2=node(teth(it,2)+1)%z-node(teth(it,4)+1)%z
    x3=(node(teth(it,3)+1)%x-node(teth(it,4)+1)%x)*1.0d2
    y3=(node(teth(it,3)+1)%y-node(teth(it,4)+1)%y)*1.1d2
    z3=node(teth(it,3)+1)%z-node(teth(it,4)+1)%z
    
    call det3(vol,x1,y1,z1,x2,y2,z2,x3,y3,z3)
    
    xp=(x-node(teth(it,4)+1)%x)*1.0d2
    yp=(y-node(teth(it,4)+1)%y)*1.1d2
    zp=z-node(teth(it,4)+1)%z

    call det3(vol1,xp,yp,zp,x2,y2,z2,x3,y3,z3)

    if( (vol*vol1 .ge. 0d0) .or. (abs(vol1/vol) .lt. 1d-12) )then

        call det3(vol2,x1,y1,z1,xp,yp,zp,x3,y3,z3)

        if( (vol*vol2 .ge. 0d0) .or. (abs(vol2/vol) .lt. 1d-12) )then

            call det3(vol3,x1,y1,z1,x2,y2,z2,xp,yp,zp)

            if( (vol*vol3 .ge. 0d0) .or. (abs(vol3/vol) .lt. 1d-12) )then

                x1=(node(teth(it,1)+1)%x-x)*1.0d2
                y1=(node(teth(it,1)+1)%y-y)*1.1d2
                z1=node(teth(it,1)+1)%z-z
                x2=(node(teth(it,2)+1)%x-x)*1.0d2
                y2=(node(teth(it,2)+1)%y-y)*1.1d2
                z2=node(teth(it,2)+1)%z-z
                x3=(node(teth(it,3)+1)%x-x)*1.0d2
                y3=(node(teth(it,3)+1)%y-y)*1.1d2
                z3=node(teth(it,3)+1)%z-z
                
                call det3(vol4,x1,y1,z1,x2,y2,z2,x3,y3,z3)

                if( (vol*vol4 .ge. 0d0) .or. (abs(vol4/vol) .lt. 1d-12) )then
                    pteth(1)%p=>node(teth(it,1)+1)
                    pteth(2)%p=>node(teth(it,2)+1)
                    pteth(3)%p=>node(teth(it,3)+1)
                    pteth(4)%p=>node(teth(it,4)+1)
                    itet=it
                    iloop=1
                    exit
                else
                    it=tethnbr(it,4)
                end if
            else
                it=tethnbr(it,3)
            end if
        else
            it=tethnbr(it,2)
        end if
    else
        it=tethnbr(it,1)
    end if
    
    if(it .eq. -1)then
        write(*,*)"Point is outside the study area! P(x,y,z):",x,y,z
        exit
    end if
end do


if(associated(pteth(1)%p) .eqv. .false.)then
    write(*,*)"point",x,y,z,"is not in tethahedron"
    stop
end if


return
end subroutine dllocat3d
!----------------------------------------------------------------------




!Update heap when the iheap element changes.s.
!-------------------------------------------------------------------------
subroutine updateheap(travelt,heap,iheap,heaptail)

implicit none
type(tstrct3d) :: travelt(maxgrid3d)
integer :: heap(maxnbnode3d)
integer :: iheap,heaptail


! For iheap which dosen't have parent node.
if(iheap/2 .eq. 0)then
    call downheap(travelt,heap,iheap,heaptail)
else if(iheap*2 .gt. heaptail)then
    call upheap(travelt,heap,iheap)
else
    if(travelt(heap(iheap/2))%t .gt. travelt(heap(iheap))%t)then
        call upheap(travelt,heap,iheap)
    else
        call downheap(travelt,heap,iheap,heaptail)
    end if
end if

return
end subroutine updateheap

!-------------------------------------------------------------------------


! Upheap
!--------------------------------------------------------------------------
subroutine upheap(travelt,heap,iheap)

implicit none
type(tstrct3d) :: travelt(maxgrid3d)
integer :: heap(maxnbnode3d)
integer :: temp
integer :: iheap
integer :: parent,child

child=iheap
parent=child/2
do while(parent .gt. 0)
    if(travelt(heap(child))%t .lt. travelt(heap(parent))%t)then
        ! Switch parent and child, include the index in %nbstat
        travelt(heap(child))%nbstat=parent
        travelt(heap(parent))%nbstat=child
        temp=heap(parent)
        heap(parent)=heap(child)
        heap(child)=temp
        child=parent
        parent=child/2
    else
        travelt(heap(child))%nbstat=child
        parent=0
    end if
end do

return
end subroutine upheap

!------------------------------------------------------------------------


! downheap
!--------------------------------------------------------------------------
subroutine downheap(travelt,heap,iheap,heaptail)

implicit none
type(tstrct3d) :: travelt(maxgrid3d)
integer :: heap(maxnbnode3d)
integer :: temp
integer :: iheap,heaptail
integer :: parent,child

parent=iheap
child=parent*2
do while(child .le. heaptail)
    
    if((child + 1) .le. heaptail)then
        if(travelt(heap(child))%t .gt. travelt(heap(child+1))%t)then
            child=child+1
        end if
    end if
    
    if(travelt(heap(parent))%t .gt. travelt(heap(child))%t)then
        ! Switch parent and child, include the index in %nbstat
        travelt(heap(child))%nbstat=parent
        travelt(heap(parent))%nbstat=child
        temp=heap(parent)
        heap(parent)=heap(child)
        heap(child)=temp
        parent=child
        child=parent*2
    else
        travelt(heap(child))%nbstat=child
        child=heaptail+1
    end if
end do

return
end subroutine downheap

!------------------------------------------------------------------------


! Determin r in which layer
!------------------------------------------------------------------------
subroutine detlayer(n,nl,layer,r)

implicit none
real(kind=8) :: layer(maxgrd1d)
real(kind=8) :: r
integer :: n,nl,i

n=0
do i=1,nl-1
    if(r .lt. layer(i))exit
    n=n+1
end do

if(r .gt. layer(nl))then
    n=n+1
end if

if((n .eq. 0) .or. (n .ge. nl))then
    write(*,*)"Error in detlayer, strct.f90!"
    write(*,*)"r is out of the layers!"
    write(*,*)"r,layer(1),layer(n): ",r,layer(1),layer(nl)
    stop
end if

return

end subroutine detlayer
!------------------------------------------------------------------------------------


! Locate p(x,y,z) in the non-uniform cubics.
subroutine locatcood3d2(n,nxl,xlayer,nyl,ylayer,nzl,zlayer,x,y,z)

implicit none
real(kind=8) :: xlayer(maxgrd1d),ylayer(maxgrd1d),zlayer(maxgrd1d)
real(kind=8) :: x,y,z
integer :: nx,ny,nz
integer :: n(8)
integer :: nxl,nyl,nzl


call detlayer(nx,nxl,xlayer,x)
call detlayer(ny,nyl,ylayer,y)
call detlayer(nz,nzl,zlayer,z)

n(1)=nxl*nyl*(nz-1)+nxl*(ny-1)+nx
n(2)=nxl*nyl*(nz-1)+nxl*(ny-1)+nx+1
n(3)=nxl*nyl*(nz-1)+nxl*ny+nx
n(4)=nxl*nyl*(nz-1)+nxl*ny+nx+1
n(5)=nxl*nyl*nz+nxl*(ny-1)+nx
n(6)=nxl*nyl*nz+nxl*(ny-1)+nx+1
n(7)=nxl*nyl*nz+nxl*ny+nx
n(8)=nxl*nyl*nz+nxl*ny+nx+1


return
end subroutine locatcood3d2
!-------------------------------------------------------------------------------


! Calculate the weighting matrix from residual vector
! npr = 1: Huber estimator
! npr = 2: Bisquare estimator
!-------------------------------------------------------------------------------
subroutine wmatrix(wei,res,nevn,npr)

implicit none
real(kind=8) :: res(maxdata3d),wei(maxdata3d)
real(kind=8) :: kwei,mar
integer :: nevn,npr,i


! Find median
if(mod(nevn,2) .eq. 0)then
        mar=(res(nevn/2)+res(nevn/2+1))/2.0d0
else
    mar=res(nevn/2+1)
end if

! npr = 1: Huber estimator
! npr = 2: Bisquare estimator
! Construct weight matrix, w(i,i)=sqrt(W(i,i))
if(npr .eq. 1)then
    kwei=1.345d0*abs(mar)/0.6745d0
    do i=1,nevn
        if(abs(res(i)) .le. kwei)then
        wei(i)=1.0d0
        else
            wei(i)=sqrt(kwei/abs(res(i)))
        end if
    end do
else if(npr .eq. 2)then
    kwei=4.685d0*abs(mar)/0.6745d0
    do i=1,nevn
        if(abs(res(i)) .le. kwei)then
            wei(i)=1.0d0-(res(i)/kwei)**2
        else
            wei(i)=0d0
        end if
    end do
else
    write(*,*)"Error!! npr is not 1 or 2!"
    stop
end if


return
end subroutine wmatrix
!-------------------------------------------------------------------------


end module strct

!-------------------------------------------------------------------------

