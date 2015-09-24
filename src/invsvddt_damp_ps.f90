program invsvddt_damp_ps


implicit none
integer,parameter :: mmax=50000,nmax=9500,nb=64
integer,parameter :: lwork=nmax+nb*nmax
real(kind=8) :: g(mmax,nmax),gtg(nmax,nmax),gtgtmp(nmax,nmax)
real(kind=8) :: d(mmax),gtd(nmax),dtmp(mmax),dtmp2(mmax),gtdtmp(nmax)
real(kind=8) :: mfinal(nmax)
real(kind=8) :: work(lwork)
real(kind=8) :: gtmp,lamdap,lamdas
real(kind=8) :: s(nmax)
real(kind=8) :: rcond
real(kind=8) :: start,finish
integer :: metagzero(nmax)
integer :: metag(nmax,2)
integer :: i,j,ir,status1,info,rank
integer :: imt,im,idmax,immax,imtmax,izero,izeromax
integer :: imtp,imts,imtpmax,imtsmax,impmax,imsmax
integer :: rowstart,tmpmetag
character(len=70) :: dfile,gfile,metagf,outf
external dgelss

call cpu_time(start)

open(31,file='invsvddt_damp_ps.inp',status='old')
read(31,*)dfile
read(31,*)metagf
read(31,*)gfile
read(31,*)outf
read(31,*)idmax,impmax,imsmax
read(31,*)lamdap,lamdas
close(31)

! Read data
open(36,file=dfile)
do i=1,idmax
    read(36,*)d(i)
end do

print *, 'Read data finished'

! Read g
immax=impmax+imsmax
metag=0
open(40,file=metagf,status='old')
rowstart=1
imt=0
imtp=0
imts=0
izero=0
do im=1,immax
    read(40,*)tmpmetag
    if(tmpmetag .ne. 0)then
        imt=imt+1
        metag(imt,1)=tmpmetag
        metag(imt,2)=rowstart
        rowstart=rowstart+metag(imt,1)
        if(im .le. impmax)then
            imtp=imtp+1
        else
            imts=imts+1
        end if
    else
        izero=izero+1
        metagzero(izero)=im
    end if
end do
imtmax=imt
imtpmax=imtp
imtsmax=imts
izeromax=izero
close(40)


open(43,file=gfile,status='old',form='unformatted',&
    &access='direct',recl=12)
g=0d0
i=0
do imt=1,imtmax
    do j=1,metag(imt,1)
        i=i+1
        read(43,rec=i,iostat=status1)gtmp,ir
        if(status1/=0)exit
        g(ir,imt)=gtmp
    end do
end do


! Normal equation and add damping
gtg=0d0
gtd=0d0
gtg(1:imtmax,1:imtmax)=matmul(transpose(g(1:idmax,1:imtmax)),&
                             &g(1:idmax,1:imtmax))
do imt=1,imtpmax
    gtg(imt,imt)=gtg(imt,imt)+lamdap
end do
do imt=imtpmax+1,imtmax
    gtg(imt,imt)=gtg(imt,imt)+lamdas
end do

gtd(1:imtmax)=matmul(transpose(g(1:idmax,1:imtmax)),d(1:idmax))
gtdtmp(1:imtmax)=gtd(1:imtmax)
gtgtmp(1:imtmax,1:imtmax)=gtg(1:imtmax,1:imtmax)

print *, 'Read g finished'
print *, 'Size of g:',idmax,imtmax

! Solve equation by svd
rcond=1d-12
call dgelss(imtmax,imtmax,1,gtg(1:imtmax,1:imtmax),imtmax,gtd(1:imtmax),&
           &imtmax,s,rcond,rank,work,lwork,info)

print *, 'Calculate model finished'
print *, 'Rank:',rank

open(51,file='sigular_values.txt',status='replace')
do im=1,immax
    write(51,*)s(im)
end do

! Out put result
mfinal=1d0
do izero=1,izeromax
    mfinal(metagzero(izero))=0d0
end do

imt=0
do im=1,immax
    if(mfinal(im) .ne. 0d0)then
        imt=imt+1
        mfinal(im)=gtd(imt)
    end if
end do

open(66,file=outf,status='replace')
do im=1,immax
    write(66,*)mfinal(im)
end do
close(66)


! Calculate residual for inversion
dtmp(1:idmax)=matmul(g(1:idmax,1:imtmax),gtd(1:imtmax))
open(75,file='invres.txt',status='replace')
do im=1,idmax
    dtmp(im)=d(im)-dtmp(im)
    write(75,*)dtmp(im)
end do
close(75)

dtmp2(1:imtmax)=matmul(gtgtmp(1:imtmax,1:imtmax),gtd(1:imtmax))
open(85,file='invres2.txt',status='replace')
do im=1,imtmax
    dtmp2(im)=gtdtmp(im)-dtmp2(im)
    write(85,*)dtmp2(im)
end do
close(85)


call cpu_time(finish)
write(*,2001)start,finish,(finish-start)/6.0d1
2001 format('start: ',f8.4,'; finish: ',f16.4,&
           &'; Time consume: ',f16.4,' min.')

stop
end

