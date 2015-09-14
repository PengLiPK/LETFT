program valmask

implicit none
real(kind=4) :: x1,y1,res,x2,y2,val
real(kind=4) :: thrshd
integer :: status1
character(len=70) :: arg1,arg2,arg3,arg4
character(len=70) :: inpf1,inpf2,outf


! Read paramters
call getarg(1,arg1)
inpf1=trim(arg1)
call getarg(2,arg2)
inpf2=trim(arg2)
call getarg(3,arg3)
outf=trim(arg3)
call getarg(4,arg4)
read(arg4,*)thrshd


! Mask val: res < thrshd
open(21,file=inpf1,status='old')
open(25,file=inpf2,status='old')
open(29,file=outf,status='replace')

do while(.true.)
    read(21,*,iostat=status1)x1,y1,res
    if(status1/=0)exit
    read(25,*)x2,y2,val
    
    if(res .gt. thrshd)then
        write(29,*)x2,y2,val
    else
        write(29,*)x2,y2,"NaN"
    end if
end do
close(21)
close(25)
close(29)

stop
end


