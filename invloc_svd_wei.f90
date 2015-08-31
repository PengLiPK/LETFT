program invloc_svd_wei

! Seperate fd matrix for veloctiy inversion and reloction by SVD

use strct
implicit none
! evnmax: max events number
! velmax: max parameters number (events paramters + vel parameters)
! idmax: max data number
! stamax: max picks per event
integer,parameter :: stamax=40,evnmax=2500,velmax=17000,idmax=70000
integer,parameter :: mmax=50,nmax=50,nb=64
integer,parameter :: lda=mmax,ldvt=nmax,ldu=mmax,lwork=nmax+nb*mmax
real(kind=8) :: a(ldu,nmax),u(ldu,nmax),vt(ldvt,nmax),s(nmax)
real(kind=8) :: g(stamax,velmax)
real(kind=8) :: tmp1(ldu),t(idmax),tly(stamax,evnmax),mt2(stamax)
real(kind=8) :: work(lwork)
real(kind=8) :: gtmp,cond,resavg,loc(4),cthrd
real(kind=8) :: res(maxdata3d),wei(maxdata3d)
integer :: metag(idmax)
integer :: nrow(evnmax),nrowtmp
integer :: ndata,nevn,nvel
integer :: i,ir,ic,id,im,ik,status1,info
integer :: icol,ily,ncol
integer :: iwork(8*mmax)
character(len=70) :: fdf,metafdf,tf,outf
external dgesdd



! Input parameters
!----------------------------------------------------

open(22,file='invloc_svd_wei.inp',status='old')
read(22,*)metafdf
read(22,*)fdf
read(22,*)tf
read(22,*)outf
read(22,*)ndata,nevn,nvel
read(22,*)cthrd
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


open(55,file=outf,status='replace')

do im=1,nevn
    if(nrow(im) .gt. 0)then

        ! Add weighting to reduce the effect of noise
        a=0d0

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

        a(1:nrow(im),1:4)=g(1:nrow(im),1:4)
        call dgesdd('All',nrow(im),4,a,lda,s,u,ldu,vt,&
                   &ldvt,work,lwork,iwork,info)
        
        ! Solve location and origin time
        tmp1(1:nrow(im))=matmul(transpose(u(1:nrow(im),1:nrow(im))),&
                               &tly(1:nrow(im),im))
        do ir=1,4
            tmp1(ir)=tmp1(ir)/s(ir)
        end do

        loc(1:4)=matmul(transpose(vt(1:4,1:4)),tmp1(1:4))

        ! Calculate residual
        mt2(1:nrow(im))=matmul(g(1:nrow(im),1:4),loc)
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



        ! Updating location and origin time
        a=0d0
        a(1:nrow(im),1:4)=g(1:nrow(im),1:4)
        call dgesdd('All',nrow(im),4,a,lda,s,u,ldu,vt,&
                   &ldvt,work,lwork,iwork,info)
        
        cond=abs(s(1)/s(4))

        ! Solve location and origin time
        tmp1(1:nrow(im))=matmul(transpose(u(1:nrow(im),1:nrow(im))),&
                               &tly(1:nrow(im),im))
        
        do ir=1,4
            tmp1(ir)=tmp1(ir)/s(ir)
        end do

        loc(1:4)=matmul(transpose(vt(1:4,1:4)),tmp1(1:4))

        ! Calculate residual for LS problem
        resavg=0d0
        do ir=5,nrow(im)
            resavg=resavg+tmp1(ir)*tmp1(ir)
        end do
        resavg=sqrt(resavg/dble(nrow(im)))

        write(*,*)im,"finished!"
    end if

    ! Output
    if(nrow(im) .gt. 0)then
        if(cond < cthrd)then
            do ir=1,4
                write(55,*)loc(ir),cond,resavg
            end do
        else
            do ir=1,4
                write(55,*)0,0,0
            end do
        end if
    else
        do ir=1,4
            write(55,*)0,0,0
        end do
    end if
end do

close(55)


stop
end
