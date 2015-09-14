! Calculate the true location
program dloc2loc

use strct
implicit none
real(kind=8) :: eloc(maxevn3d,3),dloc(maxevn3d,3) 
real(kind=8) :: dt(maxevn3d),t0
real(kind=8) :: newloc(3)
integer :: i,nevn
character(len=70) :: locfile,dlocfile,outfile




open(21,file='dloc2loc.inp',status='old')
read(21,*)locfile
read(21,*)dlocfile
read(21,*)outfile
read(21,*)nevn


open(23,file=locfile,status='old')
do i=1,nevn
    read(23,*)eloc(i,1),eloc(i,2),eloc(i,3),t0
end do
close(23)


open(26,file=dlocfile,status='old')
do i=1,nevn
    read(26,*)dt(i)
    read(26,*)dloc(i,1)
    read(26,*)dloc(i,2)
    read(26,*)dloc(i,3)
end do
close(26)


! Calculate true location
open(30,file=outfile,status='replace')
do i=1,nevn
    call loc2ear(newloc(1),newloc(2),newloc(3),&
                &eloc(i,1),eloc(i,2),eloc(i,3),&
                &-dloc(i,1),-dloc(i,2),-dloc(i,3))
    write(*,*)newloc(1),newloc(2),newloc(3),t0-dt(i)
    write(30,*)newloc(1),newloc(2),newloc(3),t0-dt(i)
end do



stop
end
