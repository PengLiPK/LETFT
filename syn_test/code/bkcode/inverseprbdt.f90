! This module defines some struct type variables which will be used in 
! main program and subroutine.
!------------------------------------------------------------------------

module typedef
    implicit none
    real(kind=8),parameter :: e=1D-12
    integer,parameter :: imodel=500000
    integer,parameter :: idata=1200000
    integer,parameter :: istanum=200

    type gstrct
        real(kind=8) :: gen
        integer :: im
        integer :: id
    end type
end module typedef
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
! This program uses conjugate gradient method solving
! least square problem.
! The equation is [G'Cd^(-1)G+lambda*Cm^(-1)]m=G'Cd^(-1)d+lambda*Cm^(-1)m0.
! Inputs are travel time, grids and frech derivative.
! Outputs are velocity models.
!--------------------------------------------------------------------------
program inverseprb

!------------------------Parameters----------------------------------------
! G:
!       "G" in [G'Cd^(-1)G+lambda*Cm^(-1)]m=G'Cd^(-1)d+lambda*Cm^(-1)m0.
! Gedge:
!       Entries of "G" coresponding to the deleted edge nodes. 
! MinLon, MaxLon, MinLan, MaxLan:
! m:
!       Model.
! d:
!       Data.
! s:
!       Parameter s in CGLS.
! r:
!       Residules.
! rabs:
!       Abusolute values of r.
! p
!       Parameter p in CGLS.
! Gp:
!        G*p.
! Alpha:
!       Parameter alpha in CGLS.
! AlphaUp:
!       Upper part of alpha.
! AlphaDown1, AlphaDown2:
!       Lower part 1 and 2 of alpha.
! Beta, Betaup, BetaDown:
!       Parameter beta in CGLS; upper and lower parts of beta.
! rmax:
!       Max value of rabs.
! m0:
!       Initial model.
! nodex, nodey:
!       Coordinates of nodes.
! edgev, edgex, edgey:
!       Velocities and coordinates of deleted edge nodes.
! lambda:
!        Damping parameter.
! cdval, cmval:
!        Covarience of data and model.
! cdinv, cminv:
!        Inverse of cdval and cmval.
! minm, maxm:
!        Maximum and minimum values of m.
! totalnodenum, edgenodenum:
!        Total number of nodes and number of deleted edge nodes.
! immax:
!        Number of nodes used in inversion.
! id, idmax:
!        Counter of data and number of data.
! im:
!        Counter of model.
! k:
!        Counter of iteration times.
! datafile:
!        Data file.
! FdFile:
!        G file.
! OputFile:
!        Output file.
! oputmodel:
!        Output model file.
! inpmodel:
!        Input model file, contains m0 and coordinates of nodes.
!--------------------------------------------------------------------------



use typedef
implicit none
type(gstrct) :: g
real(kind=8) :: m(imodel)
real(kind=8) :: d(idata)
real(kind=8) :: s(idata)
real(kind=8) :: r(imodel)
real(kind=8) :: rabs(imodel)
real(kind=8) :: p(imodel)
real(kind=8) :: Gp(idata)
real(kind=8) :: Alpha, AlphaUp, AlphaDown1,AlphaDown2
real(kind=8) :: Beta, BetaUp, BetaDown
real(kind=8) :: rmax
real(kind=8) :: nodex(imodel),nodey(imodel),nodez(imodel)
real(kind=8) :: edgev(imodel),edgex(imodel),edgey(imodel),edgez(imodel)
real(kind=8) :: lambda
real(kind=8) :: cdval,cmval
real(kind=8) :: cdinv(idata),cminv(imodel)
real(kind=8) :: minm,maxm
integer :: totalnodenum,edgenodenum
integer :: id,idmax
integer :: im,immax
integer :: ifd
integer :: k
integer :: status1
character(len=70) :: datafile
character(len=70) :: FdFile
character(len=70) :: logfile
character(len=70) :: oputmodel
character(len=70) :: inpmodel

! Input parameters
!----------------------------------------------------

open(22,file='inverseprbdt.inp')
read(22,*)FdFile
read(22,*)datafile
read(22,*)inpmodel
read(22,*)logfile
read(22,*)oputmodel
read(22,*)idmax
read(22,*)lambda
read(22,*)cdval,cmval
close(22)

!-----------------------------------------------------


! Read input model
!-----------------------------------------------------

write(*,*)cmval,cdval

open(32,file=inpmodel)
read(32,*)totalnodenum
read(32,*)edgenodenum
immax=totalnodenum-edgenodenum
if(edgenodenum .ge. 1)then
    do im=1,edgenodenum
        read(32,*)edgex(im),edgey(im),edgez(im),edgev(im)
    end do
end if

do im=1,immax
    read(32,*)nodex(im),nodey(im),nodez(im)
!    write(*,*)m0(im),im
!    imd=imd+1
!    m(im)=1.0d0/m0(im)
!    cminv(im)=1.0/cmval
end do
close(32)
m=0d0
cminv=1.0d0/cmval

!-----------------------------------------------------------


! Read data.
!------------------------------------------------------------
open(26,file=datafile)
do id=1,idmax
    read(26,*)d(id)
    cdinv(id)=1.0d0/cdval
    !write(*,*)d(id)
end do


write(*,*)maxval(cminv(1:immax)),maxval(cdinv(1:idmax))
write(*,*)" Data reading finished!"
! Reading G and Giving initial value s0, r0
!-----------------------------------------------------------

r=0d0
s=d

open(25,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=16)
ifd=0
do while(.true.)
    ifd=ifd+1
    read(25,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    s(g%id)=s(g%id)-g%gen*m(g%im)
end do

write(*,*)"s0 initial finished!"

ifd=0
do while(.true.)
    ifd=ifd+1
    read(25,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    r(g%im)=r(g%im)+g%gen*cdinv(g%id)*s(g%id)
end do

write(*,*)"r0 initial finished!"

!------------------------------------------------------------


!------------------------------------------------------------


!Initial p
!------------------------------------------------------------
!write(*,*)(cdinv(id),id=1,idmax)
!write(*,*)(cminv(im),im=1,immax)
!pause
do im=1,immax
   rabs(im)=abs(r(im))
   p(im)=r(im)
end do

rmax=maxval(rabs(1:immax))
minm=minval(m(1:immax))
maxm=maxval(m(1:immax))


write(*,*)"smax,smin: ",maxval(s(1:idmax)),minval(s(1:idmax))
write(*,*)"rmax,rmin: ",maxval(r(1:idmax)),minval(r(1:idmax))
!-------------------------------------------------------



open(27,file=logfile)
! Caculating the optimal m
!--------------------------------------------------------

k=0
    write(*,*)k
    write(27,*)k
    !write(*,*)(1/m(im),im=1,immax)
    !write(27,*)(m(im),im=1,immax)
    write(27,*)(s(im),im=1,immax)
    !write(*,*)(r(im),im=1,immax)
    !write(27,*)(r(im),im=1,immax)
    !write(27,*)(rabs(im),im=1,immax)
    !write(*,*)rmax
    write(27,*)rmax
    !pause
do while( (rmax>e) .and. (k<500)  )

    write(*,*)k
    write(27,*)k
    !write(*,*)(1/m(im),im=1,immax)
    write(27,*)(m(im),im=1,immax)
    write(27,*)(r(im),im=1,immax)
    !write(27,*)(rabs(im),im=1,immax)
    !write(*,*)rmax
    write(27,*)rmax

    k=k+1

! Caculating Alpha(k)
    Gp=0d0
    ifd=0
    do while(.true.)
        ifd=ifd+1
        read(25,rec=ifd,iostat=status1)g%gen,g%id,g%im
        if(status1/=0)exit
        Gp(g%id)=Gp(g%id)+g%gen*p(g%im)
    end do
    
    AlphaUp=0d0
    do im=1,immax
        AlphaUp=AlphaUp+r(im)*r(im)
    end do

    AlphaDown1=0d0
    do id=1,idmax
        AlphaDown1=AlphaDown1+Gp(id)*cdinv(id)*Gp(id)
    end do
    AlphaDown2=0d0
    do im=1,immax
        AlphaDown2=AlphaDown2+p(im)*cminv(im)*p(im)
    end do
    AlphaDown2=lambda*AlphaDown2
    

    Alpha=AlphaUp/(AlphaDown1+AlphaDown2)

    write(*,*)"Alpha",Alpha,AlphaUp,AlphaDown1,AlphaDown2
!    pause

! Caculating m(k+1), r(k+1), Beta(k+1), will use Gp upstairs.
    do im=1,immax
        m(im)=m(im)+Alpha*p(im)
    end do
    minm=minval(m(1:immax))
    maxm=maxval(m(1:immax))


    BetaDown=0d0
    do im=1,immax
        BetaDown=BetaDown+r(im)*r(im)
    end do


    ifd=0
    do while(.true.)
        ifd=ifd+1
        read(25,rec=ifd,iostat=status1)g%gen,g%id,g%im
        if(status1/=0)exit
        r(g%im)=r(g%im)-Alpha*g%gen*cdinv(g%id)*Gp(g%id)
    end do
    do im=1,immax
        r(im)=r(im)-Alpha*lambda*cminv(im)*p(im)
    end do


    do im=1,immax
       rabs(im)=abs(r(im))
    end do

    rmax=maxval(rabs(1:immax))

    BetaUp=0d0
    do im=1,immax
        BetaUp=BetaUp+r(im)*r(im)
    end do
    
    Beta=BetaUp/BetaDown
    write(*,*)"Beta",Beta
    write(*,*)"rmax",rmax
    write(27,*)"Beta",Beta

! Caculating p(k+1)
    do im=1,immax
        p(im)=r(im)+Beta*p(im)
    end do
end do

!-----------------------------------------------------------


close(26)
close(27)


open(29,file=oputmodel)
write(29,*)totalnodenum
write(29,*)edgenodenum
if(edgenodenum .ge. 1)then
    do im=1,edgenodenum
        write(29,*)edgex(im),edgey(im),edgez(im),edgev(im)
    end do
end if
do im=1,immax
    write(29,*)nodex(im),nodey(im),nodez(im),m(im)
end do
close(29)

! Calculate residual of Gm-d
s=d
ifd=0
do while(.true.)
    ifd=ifd+1
    read(25,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    s(g%id)=s(g%id)-g%gen*m(g%im)
end do
close(25)

open(34,file='res.txt')
do id=1,idmax
   write(34,*)s(id)
end do
close(34)

end


!----------------------------------------------------------------------------
!---------------Main program ends here---------------------------------------
!----------------------------------------------------------------------------


!############################################################################
