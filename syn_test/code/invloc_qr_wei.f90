program invloc_qr_wei

! Seperate fd matrix for veloctiy inversion and reloction by QR decomp

use strct
use solver
implicit none
integer,parameter :: evnmax=3000,velmax=5000,idmax=100000
real(kind=8) :: qt(maxrow,maxrow)
real(kind=8) :: g(maxrow,velmax,evnmax)
real(kind=8) :: t(idmax),tly(maxrow,evnmax),mt(maxrow),mt2(maxrow)
real(kind=8) :: res(maxdata3d),wei(maxdata3d)
real(kind=8) :: gtmp,gtg(4,4),gtgm(4),gtt(4)
real(kind=8) :: gr(maxrow,maxcol),loc(4)
integer :: metag(idmax)
integer :: nrow(evnmax)
integer :: ndata,nevn,nvel
integer :: i,i1,i2,i3,ir,ic,id,im,status1,iline
integer :: icol,ily,ncol
character(len=70) :: fdf,metafdf,tf,outf




! Input parameters
!----------------------------------------------------

open(22,file='invloc_qr_wei.inp',status='old')
read(22,*)metafdf
read(22,*)fdf
read(22,*)tf
read(22,*)outf
read(22,*)ndata,nevn,nvel
close(22)


! Read t
open(35,file=tf,status='old')
do id=1,ndata
    read(35,*)t(id)
end do
close(35)



! Read fd
metag=0
open(40,file=metafdf,status='old')
do id=1,ndata
    read(40,*)metag(id)
end do
close(40)



open(43,file=fdf,status='old',form='unformatted',&
    &access='direct',recl=12)
g=0d0
i=0
nrow=0
do id=1,ndata
    do im=1,metag(id)
        i=i+1
        read(43,rec=i,iostat=status1)gtmp,icol
        if(status1/=0)exit

        ! Determine the events number
        if(im .eq. 1)then
            ily=(icol/4) + 1
            nrow(ily)=nrow(ily)+1
            tly(nrow(ily),ily)=t(id)
        end if

        ! store fd matrix
        if(im .le. 4)then
            g(nrow(ily),im,ily)=gtmp
        else
            ncol=icol-nevn*4+4
            g(nrow(ily),ncol,ily)=gtmp
        end if
    end do
end do

print *, "Data reading finished!"

open(46,file=outf,status='replace')

open(50,file="gtg.txt",status='replace')
open(55,file="gtt.txt",status='replace')
i1=0
i2=0
i3=0
iline=0
do im=1,nevn
    if(nrow(im) .gt. 0)then
        ! Add weighting to lower the effects of noises
        gr=g(:,1:4,im)
        call qr(qt,gr,nrow(im),4)

        ! Calculate Q'(lower)*d and Q'(upper)*d
        mt=matmul(qt(1:nrow(im),1:nrow(im)),tly(1:nrow(im),im))

        ! Calculate dlocs and dorigin times
        loc(4)=mt(4)/gr(4,4)
        do ir=3,1,-1
            loc(ir)=mt(ir)
            do ic=4,ir+1,-1
                loc(ir)=loc(ir)-loc(ic)*gr(ir,ic)
            end do
            loc(ir)=loc(ir)/gr(ir,ir)
        end do
        
        gtg=matmul(transpose(g(1:nrow(im),1:4,im)),g(1:nrow(im),1:4,im))
        do ir=1,4
            write(50,*)gtg(ir,1),gtg(ir,2),gtg(ir,3),gtg(ir,4),im
        end do
        gtt=matmul(transpose(g(1:nrow(im),1:4,im)),tly(1:nrow(im),im))
        do ir=1,4
            write(55,*)gtt(ir),im
        end do
        gtgm=matmul(gtg(1:4,1:4),loc)

        ! Calculate residual
        mt2=matmul(g(1:nrow(im),1:4,im),loc)
        do ir=1,nrow(im)
            res(ir)=tly(ir,im)-mt2(ir)
        end do

        ! Add weighting
        call wmatrix(wei,res,nrow(im),1)
        do ir=1,nrow(im)
            do ic=1,4+nvel
                g(ir,ic,im)=g(ir,ic,im)*wei(ir)
            end do
            tly(ir,im)=tly(ir,im)*wei(ir)
        end do


        ! Update dlocations
        gr=g(:,1:4,im)
        call qr(qt,gr,nrow(im),4)

        ! Calculate Q'(lower)*d and Q'(upper)*d
        mt=matmul(qt(1:nrow(im),1:nrow(im)),tly(1:nrow(im),im))

        ! Calculate dlocs and dorigin times
        loc(4)=mt(4)/gr(4,4)
        do ir=3,1,-1
            loc(ir)=mt(ir)
            do ic=4,ir+1,-1
                loc(ir)=loc(ir)-loc(ic)*gr(ir,ic)
            end do
            loc(ir)=loc(ir)/gr(ir,ir)
        end do


        ! Calculate residual
        mt2=matmul(g(1:nrow(im),1:4,im),loc)
        do ir=1,nrow(im)
            res(ir)=tly(ir,im)-mt2(ir)
        end do

        
        write(*,*)im,"finished!"
    end if
    
    ! write output
    if(nrow(im) .gt. 0)then
        do ir=1,4
            write(46,*)loc(ir)
        end do
    else
        do ir=1,4
            write(46,*)0
        end do
    end if
end do

close(46)
close(50)
close(55)

stop
end
