! Combile s and p fd matrix and t vector
program comb_ps


implicit none

real(kind=8) :: r8
integer :: npdata,npmodel,nevn,isize
integer :: irowtmp,icol,icoltmp,status1,i,j
character(len=70) :: fdpf,tpf
character(len=70) :: fdsf,tsf
character(len=70) :: fdf,tf

open(21,file='comb_ps.inp')
read(21,*)fdpf,tpf
read(21,*)fdsf,tsf
read(21,*)npdata,npmodel,nevn
read(21,*)fdf,tf


! Combine fd
open(103,file=fdpf,status='old',form='unformatted',&
    &access='direct',recl=16)
open(105,file=fdsf,status='old',form='unformatted',&
    &access='direct',recl=16)
open(108,file=fdf,status='replace',form='unformatted',&
    &access='direct',recl=16)

i=0
do while(.true.)
    i=i+1
    read(103,rec=i,iostat=status1)r8,irowtmp,icol
    if(status1/=0)exit
    write(108,rec=i)r8,irowtmp,icol
end do
isize=i

i=0
do while(.true.)
    i=i+1
    read(105,rec=i,iostat=status1)r8,irowtmp,icol
    if(status1/=0)exit
    if(icol .gt. nevn*4)then
        icoltmp=icol+npmodel
    else
        icoltmp=icol
    end if
    j=i+isize-1
    write(108,rec=j)r8,irowtmp+npdata,icoltmp
end do

close(103)
close(105)
close(108)


! Combile t
open(113,file=tpf,status='old')
open(115,file=tsf,status='old')
open(118,file=tf,status='replace')

do while(.true.)
    read(113,*,iostat=status1)r8
    if(status1/=0)exit
    write(118,*)r8
end do

do while(.true.)
    read(115,*,iostat=status1)r8
    if(status1/=0)exit
    write(118,*)r8
end do

close(113)
close(115)
close(118)

stop
end
