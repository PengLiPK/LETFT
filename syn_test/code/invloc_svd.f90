program invloc_svd

! Seperate fd matrix for veloctiy inversion and reloction by SVD

implicit none

integer,parameter :: stamax=40,evnmax=3000,velmax=5000,idmax=100000
integer,parameter :: mmax=50,nmax=50,nb=64
integer,parameter :: lda=mmax,ldvt=nmax,ldu=mmax,lwork=nmax+nb*mmax
real(kind=8) :: a(ldu,nmax),u(ldu,nmax),vt(ldvt,nmax),s(nmax)
real(kind=8) :: g(stamax,velmax,evnmax)
real(kind=8) :: tmp1(ldu),t(idmax),tly(stamax,evnmax)
real(kind=8) :: work(lwork)
real(kind=8) :: gtmp,cond,res,loc(4),cthrd
integer :: metag(idmax)
integer :: nrow(evnmax)
integer :: ndata,nevn,nvel
integer :: i,ir,id,im,status1,info
integer :: icol,ily,ncol
integer :: iwork(8*mmax)
character(len=70) :: fdf,metafdf,tf,outf
external dgesdd



! Input parameters
!----------------------------------------------------

open(22,file='invloc_svd.inp',status='old')
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

        ! store fd matrix
        if(im .le. 4)then
            g(nrow(ily),im,ily)=gtmp
        else
            ncol=icol-nevn*4
            g(nrow(ily),ncol,ily)=gtmp
        end if
    end do
end do


open(55,file=outf,status='replace')

i=0
do im=1,nevn
    if(nrow(im) .gt. 0)then
        a=0d0
        a(1:nrow(im),1:4)=g(1:nrow(im),1:4,im)
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
        res=0d0
        do ir=5,nrow(im)
            res=res+tmp1(ir)*tmp1(ir)
        end do
        res=sqrt(res/nrow(im))

        write(*,*)im,"finished!"
    end if

    ! Output
    if(nrow(im) .gt. 0)then
        if(cond < cthrd)then
            do ir=1,4
                write(55,*)loc(ir),cond,res
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
