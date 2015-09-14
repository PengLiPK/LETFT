program invsvddt2_damp


implicit none
integer,parameter :: mmax=30000,nmax=5000,nb=64
integer,parameter :: lda=mmax,ldvt=nmax,ldu=mmax,lwork=nmax+nb*mmax
real(kind=8) :: g(mmax,nmax),gtg(nmax,nmax)
real(kind=8) :: d(mmax),tmp1(mmax),m(nmax),gtd(nmax)
real(kind=8) :: mfinal(nmax)
real(kind=8) :: work(lwork)
real(kind=8) :: rnorm,gtmp
real(kind=8) :: s(nmax),u(ldu,ldu),vt(ldvt,ldvt)
real(kind=8) :: rcond,lamda
real(kind=8) :: start,finish
integer :: metagzero(nmax)
integer :: metag(nmax,2)
integer :: i,j,ir,ic,status1,info,rank
integer :: imt,im,idmax,immax,imtmax,izero,izeromax
integer :: rowstart,tmpmetag
integer :: iwork(8*mmax)
character(len=70) :: dfile,gfile,metagf,outf
external dgesdd

call cpu_time(start)

open(31,file='invsvddt2_damp.inp',status='old')
read(31,*)dfile
read(31,*)metagf
read(31,*)gfile
read(31,*)outf
read(31,*)idmax,immax
read(31,*)lamda
close(31)

! Read data
d=0d0
open(36,file=dfile)
do i=1,idmax
    read(36,*)d(i)
end do

print *, 'Read data finished'

! Read g
metag=0
open(40,file=metagf,status='old')
rowstart=1
imt=0
izero=0
do im=1,immax
    read(40,*)tmpmetag
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


print *, 'Read g finished'
print *, 'Size of g:',idmax,imtmax


! Normal equation and add damping
gtg=0d0
gtd=0d0
gtg(1:imtmax,1:imtmax)=matmul(transpose(g(1:idmax,1:imtmax)),&
                             &g(1:idmax,1:imtmax))
do imt=1,imtmax
    gtg(imt,imt)=gtg(imt,imt)+lamda*lamda
end do

gtd(1:imtmax)=matmul(transpose(g(1:idmax,1:imtmax)),d(1:idmax))

! SVD by Lapack
call dgesdd('All',imtmax,imtmax,gtg,nmax,s,u,ldu,vt,&
           &ldvt,work,lwork,iwork,info)

print *, 'Calculate model finished'

open(51,file='sigular_values.txt',status='replace')
do im=1,immax
    write(51,*)s(im)
end do

! Calculate model
tmp1(1:imtmax)=matmul(transpose(u(1:imtmax,1:imtmax)),&
                    &gtd(1:imtmax))

do i=1,imtmax
    tmp1(i)=tmp1(i)/s(i)
end do

m(1:imtmax)=matmul(transpose(vt(1:imtmax,1:imtmax)),&
                  &tmp1(1:imtmax))


! Out put model
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

open(66,file=outf,status='replace')
do im=1,immax
    write(66,*)mfinal(im)
end do
close(66)


call cpu_time(finish)
write(*,2001)start,finish,(finish-start)/6.0d1
2001 format('start: ',f8.4,'; finish: ',f16.4,&
           &'; Time consume: ',f16.4,' min.')

stop
end

