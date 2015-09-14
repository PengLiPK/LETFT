! The program delete the tetrahedros constructed by QHULL with 0 volumn.
program del0vol


use strct
implicit none
type(srstrct3d),allocatable :: node(:)
real(kind=8) :: vol
integer :: dm
integer :: nodenum
integer :: tethnum
integer :: teth(4)
integer,allocatable :: newteth(:,:)
integer :: ind,ith
integer :: iout
character(len=70) :: nodefile
character(len=70) :: tethfile
character(len=70) :: outfile
character(len=70) :: outvolfile


! Read input file
open(21, file='del0vol.inp')
read(21,*)nodefile
read(21,*)tethfile
read(21,*)outfile
read(21,*)outvolfile

! Read the coordinates of the nodes.
open(22,file=nodefile)
read(22,*)dm
read(22,*)nodenum
allocate(node(nodenum))

do ind=1,nodenum
    read(22,*)node(ind)%x,node(ind)%y,node(ind)%z
end do
close(22)


! Read tetrahedros.
open(23,file=tethfile)

read(23,*)tethnum

open(25,file=outvolfile)

allocate(newteth(tethnum,4))
iout=0
do ith=1,tethnum
    read(23,*)teth(1),teth(2),teth(3),teth(4)

    call teth_vol(vol,node(teth(1)+1)%x,node(teth(1)+1)%y,node(teth(1)+1)%z&
                    &,node(teth(2)+1)%x,node(teth(2)+1)%y,node(teth(2)+1)%z&
                    &,node(teth(3)+1)%x,node(teth(3)+1)%y,node(teth(3)+1)%z&
                    &,node(teth(4)+1)%x,node(teth(4)+1)%y,node(teth(4)+1)%z)

    write(25,*)vol,teth(1),teth(2),teth(3),teth(4)
    if(vol .gt. 1.0e-15)then
        iout=iout+1
        newteth(iout,1)=teth(1)
        newteth(iout,2)=teth(2)
        newteth(iout,3)=teth(3)
        newteth(iout,4)=teth(4)
    end if
                    
end do

close(23)
close(25)

open(24,file=outfile)
write(24,*)iout
do ith=1,iout
    write(24,*)newteth(ith,1),newteth(ith,2),newteth(ith,3),newteth(ith,4)
end do
close(24)
stop
end
