! Add weight to G and d, used in weighted-least square problem
program addweight

implicit none
integer,parameter :: nmax=2
integer,parameter :: maxevn=200000
real(kind=8) :: f8(nmax)
real(kind=8) :: res(maxevn),absres(maxevn),wei(maxevn)
real(kind=8) :: kwei,mar,tempres
integer :: i4(nmax),npr,nevn
integer :: i,j,status1
character(len=70) :: gfile,outgf
character(len=70) :: dfile,outdf,resf


open(21,file='addweight.inp',status='old')
read(21,*)gfile
read(21,*)outgf
read(21,*)dfile
read(21,*)outdf
read(21,*)resf
read(21,*)nevn
read(21,*)npr

write(*,*)"Input file:",gfile,dfile
write(*,*)"Output file:",outgf,outdf

! Calculate k
open(30,file=resf,status='old')
do i=1,nevn
    read(30,*)res(i)
    res(i)=abs(res(i))
    absres(i)=res(i)
end do
close(30)

! sort absolute residual
do i=1,nevn-1
    do j=i,nevn
        if(res(j) .lt. res(i))then
            tempres=res(i)
            res(i)=res(j)
            res(j)=tempres
        end if
    end do
end do

! Find median
if(mod(nevn,2) .eq. 0)then
    mar=(res(nevn/2)+res(nevn/2+1))/2.0d0
else
    mar=res(nevn/2+1)
end if

! npr = 1: Huber estimator
! npr = 2: Bisquare estimator
! Construct weight matrix, w(i,i)=sqrt(W(i,i))
if(npr .eq. 1)then
    kwei=1.345d0*mar/0.6745d0
    do i=1,nevn
        if(absres(i) .le. kwei)then
            wei(i)=1.0d0
        else
            wei(i)=sqrt(kwei/absres(i))
        end if
    end do
else if(npr .eq. 2)then
    kwei=4.685d0*mar/0.6745d0
    do i=1,nevn
        if(absres(i) .le. kwei)then
            wei(i)=1.0d0-(absres(i)/kwei)**2
        else
            wei(i)=0d0
        end if
    end do
else
    write(*,*)"Error!! npr is not 1 or 2!"
end if

! Construct wG and wd
open(103,file=gfile,status='old',form='unformatted',&
    &access='direct',recl=16)
open(108,file=outgf,status='replace',form='unformatted',&
    &access='direct',recl=16)
i=0
do while(.true.)
    i=i+1
    read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2)
    if(status1/=0)exit
    f8(1)=f8(1)*wei(i4(1))
    write(108,rec=i)f8(1),i4(1),i4(2)
end do
close(103)    
close(108)

open(130,file=dfile,status='old')
open(135,file=outdf,status='replace')
do i=1,nevn
    read(130,*)f8(1)
    f8(1)=f8(1)*wei(i)
    write(135,*)f8(1)
end do
close(130)
close(135)


stop
end
