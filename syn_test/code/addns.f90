! Calculate the true location from the inversion result
program addns

implicit none
real(kind=8) :: tmpt,tmpx,tmpy,tmpz,rx,ry,rz
real(kind=8) :: sx,sy,sz
real(kind=8) :: t0,nst
integer :: i,ir,ns,nrcver,evnid
character(len=10) :: tmpstn,tmpname
character(len=70) :: inpf,outfile,nsf


open(21,file='addns.inp',status='old')
read(21,*)nsf
read(21,*)inpf
read(21,*)outfile
read(21,*)ns
close(21)


! Calculate true location
open(23,file=nsf,status='old')
open(30,file=inpf,status='old')
open(35,file=outfile,status='replace')
do i=1,ns
    read(30,*)nrcver,tmpname
    write(35,*)nrcver,tmpname
    read(30,*)tmpstn,sx,sy,sz
    write(35,2000)tmpstn,sx,sy,sz

    do ir=1,nrcver
        read(23,*)nst
        read(30,*)tmpt,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0

        write(35,2001)tmpt+nst,tmpx,tmpy,tmpz,evnid,rx,ry,rz,t0
    end do
end do
2000 format(a4,1x,f10.6,1x,f9.6,1x,f6.4)
2001 format(f10.7,1x,f10.6,1x,f9.6,1x,f6.4,1x,i4,1x,f10.6,1x,f9.6,1x,&
           &f7.4,1x,f10.7)
close(23)
close(30)
close(35)

stop
end
