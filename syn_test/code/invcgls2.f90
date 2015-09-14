! This module defines some struct type variables which will be used in 
! main program and subroutine.
!------------------------------------------------------------------------

module typedef2
    implicit none
    real(kind=8),parameter :: e=1D-12
    integer,parameter :: imodel=50000
    integer,parameter :: idata=1000000
    integer,parameter :: istanum=200

    type gstrct
        real(kind=8) :: gen
        integer :: im
        integer :: id
    end type
end module typedef2
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
! This program uses conjugate gradient method solving
! least square problem.
! The equation is [G'Cd^(-1)G+lambda*Cm^(-1)]m=G'Cd^(-1)d+lambda*Cm^(-1)m0.
! Inputs are travel time, grids and frech derivative.
! Outputs are velocity models.
!--------------------------------------------------------------------------
program invcgls2

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



use typedef2
implicit none
type(gstrct) :: g
real(kind=8) :: m(imodel),mfinal(imodel)
real(kind=8) :: d(idata)
real(kind=8) :: s(idata)
real(kind=8) :: r(imodel)
real(kind=8) :: rabs(imodel)
real(kind=8) :: p(imodel)
real(kind=8) :: Gp(idata)
real(kind=8) :: Alpha, AlphaUp, AlphaDown1,AlphaDown2
real(kind=8) :: Beta, BetaUp, BetaDown
real(kind=8) :: rmax
real(kind=8) :: lambda
real(kind=8) :: cdval,cmval
real(kind=8) :: cdinv(idata),cminv(imodel)
real(kind=8) :: minm,maxm
real(kind=8) :: thrshd,tmpdws
integer :: metag(imodel),metagzero(imodel)
integer :: imt,imtmax,izero,izeromax
integer :: id,idmax
integer :: im,immax
integer :: ifd,ifd2
integer :: k
integer :: status1
character(len=70) :: datafile
character(len=70) :: FdFile
character(len=70) :: dwsfile
character(len=70) :: oputmodel


! Input parameters
!----------------------------------------------------

open(22,file='invcgls2.inp')
read(22,*)FdFile
read(22,*)dwsfile,thrshd
read(22,*)datafile
read(22,*)oputmodel
read(22,*)idmax,immax
read(22,*)lambda
read(22,*)cdval,cmval
close(22)

!-----------------------------------------------------


! Initial model
!-----------------------------------------------------

write(*,*)cmval,cdval

m=0d0
cminv=1.0d0/cmval

!-----------------------------------------------------------


! Read data.
!------------------------------------------------------------
open(26,file=datafile,status='old')
do id=1,idmax
    read(26,*)d(id)
    cdinv(id)=1.0d0/cdval
    !write(*,*)d(id)
end do
close(26)
!-----------------------------------------------------------


! Reading DWS value by thrshd, then read fd and write new fd.
!-----------------------------------------------------------
open(30,file=dwsfile,status='old')
imt=0
izero=0
metag=0
do im=1,immax
    read(30,*)tmpdws
    if(tmpdws .ge. thrshd)then
        imt=imt+1
        metag(im)=im-izero
    else
        izero=izero+1
        metag(im)=-1
        metagzero(izero)=im
    end if
end do
imtmax=imt
izeromax=izero
close(30)

write(*,*)idmax,immax,imtmax,izeromax

open(31,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(32,file='tmpfd.dat',status='replace',form='unformatted',&
    &access='direct',recl=16)
ifd=0
ifd2=0
do while(.true.)
    ifd=ifd+1
    read(31,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    if(metag(g%im) .gt. 0)then
        ifd2=ifd2+1
        write(32,rec=ifd2)g%gen,g%id,metag(g%im)
    end if
end do
close(31)
close(32)

!-----------------------------------------------------------


write(*,*)maxval(cminv(1:imtmax)),maxval(cdinv(1:idmax))
write(*,*)" Data reading finished!"
! Reading G and Giving initial value s0, r0
!-----------------------------------------------------------

r=0d0
s=d


open(31,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(32,file='tmpfd.dat',status='replace',form='unformatted',&
    &access='direct',recl=16)
ifd=0
ifd2=0
do while(.true.)
    ifd=ifd+1
    read(31,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    if(metag(g%im) .gt. 0)then
        ifd2=ifd2+1
        write(32,rec=ifd2)g%gen,g%id,metag(g%im)
    end if
end do
close(31)
close(32)

open(35,file='tmpfd.dat',status='old',form='unformatted',&
    &access='direct',recl=16)
ifd=0
do while(.true.)
    ifd=ifd+1
    read(35,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    s(g%id)=s(g%id)-g%gen*m(g%im)
end do

write(*,*)"s0 initial finished!"

ifd=0
do while(.true.)
    ifd=ifd+1
    read(35,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    r(g%im)=r(g%im)+g%gen*cdinv(g%id)*s(g%id)
end do

write(*,*)"r0 initial finished!"

!------------------------------------------------------------


!------------------------------------------------------------


!Initial p
!------------------------------------------------------------
!pause
do im=1,imtmax
   rabs(im)=abs(r(im))
   p(im)=r(im)
end do

rmax=maxval(rabs(1:imtmax))
minm=minval(m(1:imtmax))
maxm=maxval(m(1:imtmax))


write(*,*)"smax,smin: ",maxval(s(1:idmax)),minval(s(1:idmax))
write(*,*)"rmax,rmin: ",maxval(r(1:idmax)),minval(r(1:idmax))
!-------------------------------------------------------



! Caculating the optimal m
!--------------------------------------------------------

k=0
    write(*,*)k
    !pause
do while( (rmax>e) .and. (k<500)  )

    write(*,*)k

    k=k+1

! Caculating Alpha(k)
    Gp=0d0
    ifd=0
    do while(.true.)
        ifd=ifd+1
        read(35,rec=ifd,iostat=status1)g%gen,g%id,g%im
        if(status1/=0)exit
        Gp(g%id)=Gp(g%id)+g%gen*p(g%im)
    end do
    
    AlphaUp=0d0
    do im=1,imtmax
        AlphaUp=AlphaUp+r(im)*r(im)
    end do

    AlphaDown1=0d0
    do id=1,idmax
        AlphaDown1=AlphaDown1+Gp(id)*cdinv(id)*Gp(id)
    end do
    AlphaDown2=0d0
    do im=1,imtmax
        AlphaDown2=AlphaDown2+p(im)*cminv(im)*p(im)
    end do
    AlphaDown2=lambda*AlphaDown2
    

    Alpha=AlphaUp/(AlphaDown1+AlphaDown2)

    write(*,*)"Alpha",Alpha,AlphaUp,AlphaDown1,AlphaDown2
!    pause

! Caculating m(k+1), r(k+1), Beta(k+1), will use Gp upstairs.
    do im=1,imtmax
        m(im)=m(im)+Alpha*p(im)
    end do
    minm=minval(m(1:imtmax))
    maxm=maxval(m(1:imtmax))


    BetaDown=0d0
    do im=1,imtmax
        BetaDown=BetaDown+r(im)*r(im)
    end do


    ifd=0
    do while(.true.)
        ifd=ifd+1
        read(35,rec=ifd,iostat=status1)g%gen,g%id,g%im
        if(status1/=0)exit
        r(g%im)=r(g%im)-Alpha*g%gen*cdinv(g%id)*Gp(g%id)
    end do
    do im=1,imtmax
        r(im)=r(im)-Alpha*lambda*cminv(im)*p(im)
    end do


    do im=1,imtmax
       rabs(im)=abs(r(im))
    end do

    rmax=maxval(rabs(1:imtmax))

    BetaUp=0d0
    do im=1,imtmax
        BetaUp=BetaUp+r(im)*r(im)
    end do
    
    Beta=BetaUp/BetaDown
    write(*,*)"Beta",Beta
    write(*,*)"rmax",rmax

! Caculating p(k+1)
    do im=1,imtmax
        p(im)=r(im)+Beta*p(im)
    end do
end do

!-----------------------------------------------------------


! Write output file
!-----------------------------------------------------------
mfinal=1d0
do izero=1,izeromax
    mfinal(metagzero(izero))=0d0
end do

imt=0
do im=1,immax
    if(mfinal(im) .ne. 0d0)then
        imt=imt+1
        mfinal(im)=m(imt)
    end if
end do

open(49,file=oputmodel)
do im=1,immax
    write(49,*)mfinal(im)
end do
close(49)
!-----------------------------------------------------------


! Calculate residual of Gm-d
s=d
ifd=0
do while(.true.)
    ifd=ifd+1
    read(35,rec=ifd,iostat=status1)g%gen,g%id,g%im
    if(status1/=0)exit
    s(g%id)=s(g%id)-g%gen*m(g%im)
end do
close(35)

open(54,file='res.txt')
do id=1,idmax
   write(54,*)s(id)
end do
close(54)

end


!----------------------------------------------------------------------------
!---------------Main program ends here---------------------------------------
!----------------------------------------------------------------------------


!############################################################################
