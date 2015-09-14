program fdsep_svd

! Seperate fd matrix for veloctiy inversion and reloction by SVD

implicit none

integer,parameter :: stamax=40,evnmax=3000,velmax=5000,idmax=100000
integer,parameter :: mmax=50,nmax=50,nb=64
integer,parameter :: lda=mmax,ldvt=nmax,ldu=mmax,lwork=nmax+nb*mmax
real(kind=8) :: u(ldu,nmax),vt(ldvt,nmax)
real(kind=8) :: g(stamax,velmax,evnmax),mg(stamax,velmax)
real(kind=8) :: t(idmax),tly(stamax,evnmax),mt(stamax)
real(kind=8) :: work(lwork)
real(kind=8) :: s(nmax)
real(kind=8) :: gtmp
integer :: metag(idmax)
integer :: nrow(evnmax)
integer :: ndata,nevn,nvel
integer :: i,ir,ic,id,im,status1,info
integer :: icol,ily,ncol
integer :: iwork(8*mmax)
character(len=70) :: fdf,metafdf,tf,outfdf,outtf
external dgesvd



! Input parameters
!----------------------------------------------------

open(22,file='fdsep_svd.inp',status='old')
read(22,*)metafdf
read(22,*)fdf
read(22,*)tf
read(22,*)outfdf
read(22,*)outtf
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
            ncol=icol-nevn*4
            g(nrow(ily),ncol,ily)=gtmp
        end if
    end do
end do


open(50,file=outfdf,status='replace',form='unformatted',&
     &access='direct',recl=16)
open(55,file=outtf,status='replace')

i=0
do im=1,nevn
    if(nrow(im) .gt. 0)then
        call dgesdd('A',nrow(im),4,g(:,1:4,im),lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)
        !do ir=1,nrow(im)
        !    do ic=1,nrow(im)
        !        ut(ir,ic)=u(ic,ir)
        !    end do
        !end do
        
        do ir=1,nrow(im)-4
            do ic=1,nvel
                mg(ir,ic)=0d0
                do id=1,nrow(im)
                    mg(ir,ic)=mg(ir,ic)+u(id,ir+4)*g(id,ic+4,im)
                end do
                i=i+1
                write(50,rec=i)mg(ir,ic),ir,ic
            end do
        end do

        do ir=1,nrow(im)-4
            mt(ir)=0d0
            do id=1,nrow(im)
                mt(ir)=mt(ir)+u(id,ir+4)*tly(id,im)
            end do
            write(55,*)mt(ir)
        end do
        write(*,*)im,"finished!"
    end if
end do

close(50)
close(55)


stop
end
