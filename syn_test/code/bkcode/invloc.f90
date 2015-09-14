! Compute locations and origin times
program invloc

implicit none
integer,parameter :: maxrow=4,maxcol=2000,maxdata=50000
real(kind=8) :: r(maxrow,maxrow),qa(maxrow,maxcol)
real(kind=8) :: dt(maxdata),dsl(maxcol),locall(maxdata)
real(kind=8) :: loc(4),dt2(4),dt3(4)
real(kind=8) :: rtmp,qatmp
integer :: evnid(maxdata),evnpara(maxdata)
integer :: ndata,nevnall,nsl,nevn
integer :: ievn,irow,icol,i,j,i1,i2,ir1,id,status1
character(len=70) :: fdfile,dtfile,rfile,dslfile,outfile


open(21,file='invloc.inp')
read(21,*)fdfile
read(21,*)dtfile
read(21,*)rfile
read(21,*)dslfile
read(21,*)outfile
read(21,*)nevnall,ndata,nsl


! Read dt
open(101,file=dtfile,status='old')
do id=1,ndata
    read(101,*)dt(id)
end do
close(101)

! Read dsl
open(102,file=dslfile,status='old')
do id=1,nsl
    read(102,*)dsl(id)
end do
close(102)


! Calculate the number of events
if(mod(ndata,4) .ne. 0)then
    write(*,*)"Error!! Data number can't divided by 4.",ndata
    stop
else
    nevn=ndata/4
end if


open(103,file=rfile,status='old',form='unformatted',&
    &access='direct',recl=20)
open(105,file=fdfile,status='old',form='unformatted',&
    &access='direct',recl=16)

i1=0
i2=0
do id=1,nevn

    
    ! Read R matrix, R is 4x4
    ir1=1
    r=0d0
    do while(.true.)
        i1=i1+1
        read(103,rec=i1,iostat=status1)rtmp,irow,icol,ievn
        if(status1 .ne. 0 .and. id .ne. nevn)then
            print *, "nevn /= evn number from R-matrix file!"
            exit
        else if(status1 .ne. 0)then
            exit
        else
            if(irow .lt. ir1)then
                i1=i1-1
                exit
            else
                ir1=irow
                evnid(id)=ievn
                r(irow,icol)=rtmp
            end if
        end if
    end do

    
    ! Read Q'(upper)*A2
    ir1=1
    qa=0d0
    do while(.true.)
        i2=i2+1
        read(105,rec=i2,iostat=status1)qatmp,irow,icol
        if(status1 .ne. 0 .and. id .ne. nevn)then
            print *, "nevn /= evn number from fd file!"
            exit
        else if(status1 .ne. 0)then
            exit
        else
            if(irow .lt. ir1)then
                i2=i2-1
                exit
            else
                ir1=irow
                qa(irow,icol)=qatmp
            end if
        end if
    end do

    ! Calculate qa*dsl
    dt2=matmul(qa(1:4,1:nsl),dsl(1:nsl))

    ! Calculate dlocs and dorigin times
    do i=1,4
        dt3(i)=dt(i+(id-1)*4)-dt2(i)
    end do

    loc(4)=dt3(4)/r(4,4)
    do i=3,1,-1
        loc(i)=dt3(i)
        do j=4,i+1,-1
            loc(i)=loc(i)-loc(j)*r(i,j)
        end do
        loc(i)=loc(i)/r(i,i)
    end do

    do i=1,4
        locall((id-1)*4+i)=loc(i)
    end do

end do

close(103)
close(105)

! Write outputfile
evnpara=0
do id=1,nevn
    evnpara(evnid(id))=1
end do

open(120,file=outfile,status='replace')
i1=0
do id=1,nevnall
    if(evnpara(id) .eq. 0)then
        do i=1,4
            write(120,*)0d0
        end do
    else
        do i=1,4
            i1=i1+1
            write(120,*)locall(i1)
        end do
    end if
end do
close(120)

stop
end
