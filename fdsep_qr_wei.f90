program fdsep_qr_wei

! Seperate fd matrix for veloctiy inversion and reloction by QR decomp

use solver
use strct
implicit none
! evnmax: max events number
! velmax: max parameters number (events paramters + 4)
! idmax: max data number
integer,parameter :: evnmax=2500,velmax=18000,idmax=70000
real(kind=8) :: qt(maxrow,maxrow)
real(kind=8) :: g(maxrow,velmax),mg(maxrow,velmax)
real(kind=8) :: t(idmax),tly(maxrow,evnmax),mt(maxrow),mt2(maxrow)
real(kind=8) :: res(maxdata3d),wei(maxdata3d)
real(kind=8) :: gtmp
real(kind=8) :: gr(maxrow,maxcol),loc(4)
integer :: metag(idmax)
integer :: nrow(evnmax),nrowtmp
integer :: ndata,nevn,nvel
integer :: i,i1,i2,i3,ir,ic,id,im,ik,status1,iline
integer :: icol,ily,ncol
character(len=70) :: fdf,metafdf,tf,outfdf1,outtf1,outrf,outfdf2,outtf2




! Input parameters
!----------------------------------------------------

open(22,file='fdsep_qr_wei.inp',status='old')
read(22,*)metafdf
read(22,*)fdf
read(22,*)tf
read(22,*)outfdf1
read(22,*)outtf1
read(22,*)outrf
read(22,*)outfdf2
read(22,*)outtf2
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

    end do
end do

print *, "Data reading finished!"

open(45,file=outfdf1,status='replace',form='unformatted',&
     &access='direct',recl=16)
open(46,file=outtf1,status='replace')
open(47,file=outrf,status='replace',form='unformatted',&
     &access='direct',recl=20)
open(50,file=outfdf2,status='replace',form='unformatted',&
     &access='direct',recl=16)
open(55,file=outtf2,status='replace')

i1=0
i2=0
i3=0
iline=0
do im=1,nevn
    if(nrow(im) .gt. maxrow)then
        write(*,*)nrow(im),maxrow
        write(*,*)"Parameter maxrow in solver.f90 is too small,&
                 &replace it with a larger one."
        stop
    else if(nrow(im) .gt. 0)then

        ! Add weighting to lower the effects of noises
        i=0
        nrowtmp=0
        g=0d0
        do id=1,ndata
            do ik=1,metag(id)
                i=i+1
                read(43,rec=i,iostat=status1)gtmp,icol
                if(status1/=0)exit
                if(ik .eq. 1)then
                    ily=(icol/4)+1
                end if
                
                if(ily .eq. im)then
                    if(ik .eq. 1)then
                        nrowtmp=nrowtmp+1
                    end if
                    if(ik .le. 4)then
                        g(nrowtmp,ik)=gtmp
                    else
                        ncol=icol-nevn*4+4
                        g(nrowtmp,ncol)=gtmp
                    end if
                end if
            end do
        end do
        
        gr=g(:,1:4)
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
        mt2=matmul(g(1:nrow(im),1:4),loc)
        do ir=1,nrow(im)
            res(ir)=tly(ir,im)-mt2(ir)
        end do

        ! Add weighting
        call wmatrix(wei,res,nrow(im),1)
        do ir=1,nrow(im)
            do ic=1,4+nvel
                g(ir,ic)=g(ir,ic)*wei(ir)
            end do
            tly(ir,im)=tly(ir,im)*wei(ir)
        end do



        ! Parameter seperation
        gr=g(:,1:4)
        call qr(qt,gr,nrow(im),4)
        
        ! Calculate Q'(lower)*A and Q'(upper)*A
        mg=matmul(qt(1:nrow(im),1:nrow(im)),g(1:nrow(im),5:4+nvel))
        do ir=1,nrow(im)
            if(ir .gt. 4)then
                ! Store Q'(lower)*A for slowness inversion
                iline=iline+1
                do ic=1,nvel
                    if(mg(ir,ic) .ne. 0d0)then
                        i1=i1+1
                        write(50,rec=i1)mg(ir,ic),iline,ic
                    end if
                end do
            else
                ! Store Q'(upper)*A for relocation
                do ic=1,nvel
                    if(mg(ir,ic) .ne. 0d0)then
                        i2=i2+1
                        write(45,rec=i2)mg(ir,ic),ir,ic
                    end if
                end do
            end if
        end do

        ! Store non-zero elements of R for relocation
        do ir=1,4
            do ic=1,4
                if(gr(ir,ic) .ne. 0d0)then
                    i3=i3+1
                    write(47,rec=i3)gr(ir,ic),ir,ic,im
                end if
            end do
        end do

        ! Calculate Q'(lower)*d and Q'(upper)*d
        mt=matmul(qt(1:nrow(im),1:nrow(im)),tly(1:nrow(im),im))
        do ir=1,nrow(im)
            if(ir .gt. 4)then
                write(55,*)mt(ir)
            else
                write(46,*)mt(ir)
            end if
        end do
        write(*,*)im,"finished!"
    end if
end do

close(45)
close(46)
close(47)
close(50)
close(55)


stop
end
