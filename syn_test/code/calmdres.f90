program calmdres


implicit none
integer,parameter :: mmax=27000,nmax=5000,nb=64
integer,parameter :: lda=mmax,ldvt=nmax,ldu=mmax,lwork=nmax+nb*mmax
real(kind=8) :: g(lda,nmax),abk(lda,nmax),u(ldu,nmax),vt(ldvt,nmax)
real(kind=8) :: ainv(nmax,lda)
real(kind=8) :: mfinal(nmax),r(nmax)
real(kind=8) :: dwsfinal(nmax),dws(nmax)
real(kind=8) :: work(lwork)
real(kind=8) :: s(nmax)
real(kind=8) :: gtmp
integer :: rctfinal(nmax),rct(nmax)
integer :: metagzero(nmax)
integer :: metag(nmax,2)
integer :: i,j,k,ir,status1,info
integer :: imt,im,idmax,immax,imtmax,izero,izeromax
integer :: rowstart,tmpmetag
integer :: iwork(8*mmax)
character(len=70) :: gfile,metagf,outrf,outrctf,outdwsf
!external dnrm2
external dgesdd


open(31,file='calmdres.inp',status='old')
read(31,*)metagf
read(31,*)gfile
read(31,*)outrf
read(31,*)outrctf
read(31,*)outdwsf
read(31,*)idmax,immax
close(31)


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
dws=0d0
rct=0
i=0
do imt=1,imtmax
    do j=1,metag(imt,1)
        i=i+1
        read(43,rec=i,iostat=status1)gtmp,ir
        if(status1/=0)exit
        g(ir,imt)=gtmp
        dws(imt)=dws(imt)+gtmp
        if(gtmp .gt. 0)then
            rct(imt)=rct(imt)+1
        end if
    end do
end do


abk=g
print *, 'Read g finished'
print *, 'immax=',immax,',  imtmax=',imtmax

call dgesvd('S','A',idmax,imtmax,g,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)

print *, 'Calculate model finished'

! Calculate model resolution matrix
do i=1,imtmax
    do j=1,imtmax
        vt(i,j)=vt(i,j)/s(i)
    end do
end do

do i=1,imtmax
    do j=1,idmax
        ainv(i,j)=0d0
        do k=1,imtmax
            ainv(i,j)=ainv(i,j)+vt(k,i)*u(j,k)
        end do
   end do
end do

! Only calculate the diagnal entries
do i=1,imtmax
    r(i)=0d0
    do k=1,idmax
        r(i)=r(i)+ainv(i,k)*abk(k,i)
    end do
end do


open(51,file='sigular_values.txt',status='replace')
do im=1,immax
    write(51,*)s(im)
end do


! Calculate final r,rct,dws
mfinal=1d0
rctfinal=1
dwsfinal=1d0
do izero=1,izeromax
    mfinal(metagzero(izero))=0d0
    rctfinal(metagzero(izero))=0
    dwsfinal(metagzero(izero))=0d0
end do

imt=0
do im=1,immax
    if(mfinal(im) .ne. 0d0)then
        imt=imt+1
        mfinal(im)=r(imt)
        rctfinal(im)=rct(imt)
        dwsfinal(im)=dws(imt)
    end if
end do

! Write output files
open(66,file=outrf,status='replace')
do im=1,immax
    write(66,*)mfinal(im)
end do
close(66)

open(70,file=outrctf,status='replace')
do im=1,immax
    write(70,*)rctfinal(im)
end do
close(70)

open(75,file=outdwsf,status='replace')
do im=1,immax
    write(75,*)dwsfinal(im)
end do
close(75)


stop
end

