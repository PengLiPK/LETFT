! This module defines some struct type variables which will be used in 
! main program and subroutine.
!------------------------------------------------------------------------

module typedef
    implicit none
    real(kind=8),parameter :: e=1D-12
    integer,parameter :: imodel=400000
    integer,parameter :: idata=1200000
    integer,parameter :: istanum=200

    type gstrct
        real(kind=8) :: gen
        integer :: row
    end type

contains

! Calculate 2-norm of a vector.
!------------------------------------------------------------------------
subroutine norm2(vnorm,vector,length)

implicit none
real(kind=8) :: vnorm
real(kind=8) :: vector(:)
integer :: length,i


vnorm=0d0
do i=1,length
    vnorm=vnorm+vector(i)*vector(i)
end do
vnorm=sqrt(vnorm)

return
end subroutine norm2


! Calculate the product of two vectors in sparse matrix.
!--------------------------------------------------------------------------
subroutine vectormult(vpdt,v1,nv1,v2,nv2)

implicit none
type(gstrct) :: v1(:),v2(:)
real(kind=8) :: vpdt
integer :: nv1,nv2
integer :: i,j,k


vpdt=0d0

k=1
do i=1,nv1
    do j=k,nv2
        if(v1(i)%row .lt. v2(j)%row)then
            exit
        else if(v1(i)%row .eq. v2(j)%row)then
            vpdt=vpdt+v1(i)%gen*v2(j)%gen
            k=j+1
            exit
        end if
    end do
    if(k .gt. nv2)exit
end do

return
end subroutine vectormult
!------------------------------------------------------------------------


! Symmetric Gauss-Seidel preconditioner solver.
! Preconditioner M=(D+L)D^-1(D+U)
! Solver:
! (D+L)*p=r
! D^-1(D+U)*p=p
!-------------------------------------------------------------------------
subroutine sgssol(p,g,metag,r,immax)

implicit none
type(gstrct) :: g(:)
real(kind=8) :: r(imodel)
real(kind=8) :: p(imodel)
real(kind=8) :: mdiag(imodel)
real(kind=8) :: m1,m2
integer :: metag(imodel,2)
integer :: immax
integer :: i,j,m


! Solve (D+L)*p=r
p=r
m1=0d0
do i=1,metag(1,1)
    m1=m1+g(i)%gen*g(i)%gen
end do
p(1)=r(1)/m1

mdiag(1)=m1

do i=2,immax
    do j=1,i-1
        ! Calculate m1(i,j)
        call vectormult(m1,g(metag(i,2):metag(i,2)+metag(i,1)-1),metag(i,1),&
                       &g(metag(j,2):metag(j,2)+metag(j,1)-1),metag(j,1))
        p(i)=p(i)-m1*p(j)
    end do

    m1=0d0
    do m=1,metag(i,1)
        m1=m1+g(m+metag(i,2)-1)%gen*g(m+metag(i,2)-1)%gen
    end do
    mdiag(i)=m1
    p(i)=p(i)/m1
end do

! Solve D^-1(D+U)*p=p
do i=immax-1,1,-1
    do j=i+1,immax
        call vectormult(m2,g(metag(i,2):metag(i,2)+metag(i,1)-1),metag(i,1),&
                       &g(metag(j,2):metag(j,2)+metag(j,1)-1),metag(j,1))
        m2=m2/mdiag(i)
        p(i)=p(i)-m2*p(j)
    end do
    p(i)=p(i)
end do

return
end subroutine sgssol
!-------------------------------------------------------------------------


end module typedef
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
! This program uses conjugate gradient method solving
! least square problem.
! The equation is [G'Cd^(-1)G+lambda*Cm^(-1)]m=G'Cd^(-1)d+lambda*Cm^(-1)m0.
! Inputs are travel time, grids and frech derivative.
! Outputs are velocity models.
!--------------------------------------------------------------------------
program inversepcg

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
type(gstrct),allocatable :: g(:)
real(kind=8) :: m(imodel)
real(kind=8) :: mfinal(imodel)
real(kind=8) :: d(idata)
real(kind=8) :: s(idata)
real(kind=8) :: r(imodel)
real(kind=8) :: p(imodel)
real(kind=8) :: z(imodel)
real(kind=8) :: gtg,gtgp(imodel)
real(kind=8) :: Alpha, AlphaUp, AlphaDown1,AlphaDown2
real(kind=8) :: Beta, BetaUp, BetaDown
real(kind=8) :: rmax
real(kind=8) :: m0(imodel),nodex(imodel),nodey(imodel),nodez(imodel)
real(kind=8) :: edgev(imodel),edgex(imodel),edgey(imodel),edgez(imodel)
real(kind=8) :: lambda
real(kind=8) :: cdval,cmval
real(kind=8) :: cdinv(idata),cminv(imodel)
real(kind=8) :: minm,maxm
real(kind=8) :: start,finish
real(kind=8) :: rtime(8)
integer :: tmpmetag
integer :: metagzero(imodel)
integer :: metag(imodel,2)
integer :: rowstart
integer :: totalnodenum,edgenodenum
integer :: id,idmax
integer :: im,immax
integer :: imt,imtmax
integer :: izero,izeromax
integer :: ifd,ia
integer :: k
character(len=70) :: datafile
character(len=70) :: FdFile
character(len=70) :: fdmetafile
character(len=70) :: logfile
character(len=70) :: oputmodel
character(len=70) :: inpmodel


call cpu_time(start)

! Input parameters
!----------------------------------------------------

open(22,file='invpcgdt.inp')
read(22,*)FdFile
read(22,*)fdmetafile
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

open(23,file=inpmodel)
read(23,*)totalnodenum
read(23,*)edgenodenum
immax=totalnodenum-edgenodenum
if(edgenodenum .ge. 1)then
    do im=1,edgenodenum
        read(23,*)edgex(im),edgey(im),edgez(im),edgev(im)
    end do
end if

do im=1,immax
    read(23,*)nodex(im),nodey(im),nodez(im),m0(im)
!    m(im)=1.0d0/m0(im)
end do
m=0d0
close(23)
cminv=1.0d0/cmval

!-----------------------------------------------------------


! Read data.
!------------------------------------------------------------
open(26,file=datafile)
do id=1,idmax
    read(26,*)d(id)
    cdinv(id)=1.0d0/cdval
end do
close(26)

write(*,*)maxval(cminv(1:immax)),maxval(cdinv(1:idmax))
write(*,*)" Data reading finished!"

! Read metag and G
!------------------------------------------------------------
metag=0
open(42,file=fdmetafile,status='old')
rowstart=1
imt=0
izero=0
do im=1,immax
    read(42,*)tmpmetag
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
close(42)

write(*,*)rowstart-1,imtmax,izeromax

allocate(g(rowstart-1))

open(45,file=FdFile,status='old',form='unformatted',&
    &access='direct',recl=12)

do ifd=1,rowstart-1
    read(45,rec=ifd)g(ifd)%gen,g(ifd)%row
end do


! Giving initial value s0, r0
!-----------------------------------------------------------

r=0d0

do im=1,imtmax
    do ifd=1,metag(im,1)
        r(im)=r(im)+g(ifd+metag(im,2)-1)%gen*&
             &cdinv(g(ifd+metag(im,2)-1)%row)*&
             &d(g(ifd+metag(im,2)-1)%row)
    end do
end do

write(*,*)"r0 initial finished!"
write(*,*)"rmax,rmin: ",maxval(r(1:idmax)),minval(r(1:idmax))

!------------------------------------------------------------



!Initial p
!------------------------------------------------------------
call cpu_time(rtime(1))
call sgssol(p,g,metag,r,imtmax)
call cpu_time(rtime(2))
write(*,*)rtime(2)-rtime(1),'seconds'

z=p;

minm=minval(m(1:imtmax))
maxm=maxval(m(1:imtmax))


write(*,*)"smax,smin: ",maxval(s(1:idmax)),minval(s(1:idmax))
!-------------------------------------------------------



open(55,file=logfile)
! Caculating the optimal m
!--------------------------------------------------------
k=0
call norm2(rmax,r,imtmax)
do while( (rmax>e) .and. (k<300)  )

    write(*,*)k
    write(55,*)k
    write(55,*)rmax
    write(*,*)"Norm2 of res:",rmax

    k=k+1

! Caculating Alpha(k)
    
    AlphaUp=0d0
    do im=1,imtmax
        AlphaUp=AlphaUp+r(im)*z(im)
    end do

    AlphaDown1=0d0
    gtgp=0d0
    do im=1,imtmax
        do ia=1,imtmax
            ! Calculate gtg(im,ia)
            call vectormult(gtg,g(metag(im,2):metag(im,2)+metag(im,1)-1),metag(im,1),&
                           &g(metag(ia,2):metag(ia,2)+metag(ia,1)-1),metag(ia,1))
            gtgp(im)=gtgp(im)+gtg*p(ia)
        end do
    end do

    do im=1,imtmax
        AlphaDown1=AlphaDown1+p(im)*gtgp(im)
    end do
    AlphaDown2=0d0
    do im=1,imtmax
        AlphaDown2=AlphaDown2+p(im)*cminv(im)*p(im)
    end do
    AlphaDown2=lambda*AlphaDown2
    

    Alpha=AlphaUp/(AlphaDown1+AlphaDown2)

    write(*,*)"Alpha",Alpha,AlphaUp,AlphaDown1,AlphaDown2

! Update m(k+1), r(k+1), z(k+1), Beta(k+1), will use gtgp upstairs.
    do im=1,imtmax
        m(im)=m(im)+Alpha*p(im)
    end do
    minm=minval(m(1:imtmax))
    maxm=maxval(m(1:imtmax))


    BetaDown=0d0
    do im=1,imtmax
        BetaDown=BetaDown+r(im)*z(im)
    end do


    do im=1,imtmax
        r(im)=r(im)-Alpha*gtgp(im)
    end do
    do im=1,imtmax
        r(im)=r(im)-Alpha*lambda*cminv(im)*p(im)
    end do

    call norm2(rmax,r,imtmax)

    call sgssol(z,g,metag,r,imtmax)

    BetaUp=0d0
    do im=1,imtmax
        BetaUp=BetaUp+r(im)*z(im)
    end do
    
    Beta=BetaUp/BetaDown
    write(*,*)"Beta",Beta
    write(55,*)"Beta",Beta

! Update p(k+1)
    do im=1,imtmax
        p(im)=z(im)+Beta*p(im)
    end do

end do

close(55)
!--------------------------------------------------------------------------


! Output velocity model
!----------------------------------------------------------------------------

mfinal=1d0
do izero=1,izeromax 
    mfinal(metagzero(izero))=0d0
end do

imt=0
do im=1,immax
    if(mfinal(im) .ne. 0d0)then
        imt=imt+1
        mfinal(im)=1.0d0/m(imt)
    end if
end do


open(66,file=oputmodel)
write(66,*)totalnodenum
write(66,*)edgenodenum
if(edgenodenum .ge. 1)then
    do im=1,edgenodenum
        write(66,*)edgex(im),edgey(im),edgez(im),0
    end do
end if
do im=1,immax
    write(66,*)nodex(im),nodey(im),nodez(im),mfinal(im)
end do
close(66)
!-----------------------------------------------------------------------------


! Calculate and output residual of Gm-d
! Gm=G(:,1)*m1+G(:,2)*m2+...+G(:,k)*mk
!--------------------------------------------------------------------------------
do im=1,imtmax
    do ifd=1,metag(im,1)
        d(g(ifd+metag(im,2)-1)%row)=d(g(ifd+metag(im,2)-1)%row)-&
                                   &g(ifd+metag(im,2)-1)%gen*m(im)
    end do
end do

call norm2(rmax,d,idmax)

write(*,*)"||Gm-d||=",rmax

open(75,file='res.txt')
do id=1,idmax
   write(75,*)d(id)
end do
write(75,*)rmax
close(75)
!-------------------------------------------------------------------------------


call cpu_time(finish)
write(*,2000)start,finish,finish-start
2000 format('Start: ',f8.4,'; Finish: ',f16.4,&
           &'; Time consume: ',f16.4,'sec.')

stop

end


!--------------------------------------------------------------------------------
!---------------Main program ends here-------------------------------------------
!--------------------------------------------------------------------------------


!############################################################################
