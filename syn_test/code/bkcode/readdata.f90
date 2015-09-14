program readdata
!#################################################################################################
! Read binary files in adaptive tomography program.
! ftype = 1: Output of sourceXXXXX.out
!       = 2: Ray path file
!       = 3: fd file
!       = 4: fdcsm file
!       = 5: r file
!#################################################################################################

implicit none
integer,parameter :: nmax=100 
real(kind=8) :: f8(nmax)
real(kind=4) :: f4(nmax)
integer :: i4(nmax)
integer :: ftype
integer :: i,j,status1
character(len=70) :: dfile,outf


open(21,file='readdata.inp',status='old')
read(21,*)dfile
read(21,*)outf
read(21,*)ftype

write(*,*)"Input file:",dfile
write(*,*)"Output file:",outf
write(*,*)"file type:",ftype

if(ftype .eq. 1)then
    open(103,file=dfile,status='old',form='unformatted',&
        &access='direct',recl=60)
    open(108,file=outf,status='replace')
    i=0
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),f8(2),f8(3),f8(4),f8(5),f8(6),i4(1),i4(2),i4(3)
        if(status1/=0)exit
        write(108,1001)f8(1),f8(2),f8(3),f8(4),f8(5),f8(6),i4(1),i4(2),i4(3)
1001 format(1x,f16.12,1x,f16.12,1x,f16.12,1x,f16.12,1x,f16.12,1x,f16.12,1x,i6.6,1x,i5.5,1x,i1)
    end do
    close(103)    
    close(108)
else if(ftype .eq. 2)then
    open(103,file=dfile,status='old',form='unformatted',&
        &access='direct',recl=24)
    open(108,file=outf,status='replace')
    i=0
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)i4(1)
        if(status1/=0)exit
        write(108,1002)i4(1)
1002 format('N ',i6)
        do j=1,i4(1)
            i=i+1
            read(103,rec=i,iostat=status1)f8(1),f8(2),f8(3)
            write(108,1003)f8(1),f8(2),f8(3)
1003 format(1x,f16.12,1x,f16.12,1x,f16.12)
        end do
    end do
    close(103)    
    close(108)
else if(ftype .eq. 3)then
    open(103,file=dfile,status='old',form='unformatted',&
        &access='direct',recl=16)
    open(108,file=outf,status='replace')
    i=0
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2)
        if(status1/=0)exit
        write(108,1004)f8(1),i4(1),i4(2)
1004 format(f16.12,1x,i6,1x,i5)
    end do
    close(103)    
    close(108)
else if(ftype .eq. 4)then
    open(103,file=dfile,status='old',form='unformatted',&
        &access='direct',recl=12)
    open(108,file=outf,status='replace')
    i=0
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),i4(1)
        if(status1/=0)exit
        write(108,1005)f8(1),i4(1)
1005 format(f16.12,1x,i6)
    end do
    close(103)    
    close(108)
else if(ftype .eq. 5)then
    open(103,file=dfile,status='old',form='unformatted',&
        &access='direct',recl=20)
    open(108,file=outf,status='replace')
    i=0
    do while(.true.)
        i=i+1
        read(103,rec=i,iostat=status1)f8(1),i4(1),i4(2),i4(3)
        if(status1/=0)exit
        write(108,1006)f8(1),i4(1),i4(2),i4(3)
1006 format(f16.12,1x,i6,1x,i6,1x,i6)
    end do
    close(103)    
    close(108)
end if


stop
end
