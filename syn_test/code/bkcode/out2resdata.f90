program out2resdata


implicit none



open(21,file='out2resdata.inp',status='old')
read(21,*)resulf
read(21,*)dataf
read(21,*)ievn,ivel


open(22,file=resulf,status='old')
do i=1,ievn
    read(22,*)dt(i)
    read(22,*)dx(i)
    read(22,*)dy(i)
    read(22,*)dz(i)
end do
do i=1,ivel
    read(22,*)vel(i)
end do
close(22)


open(23,file=dataf,status='old')
ndata=0
ttlevn=0
do is=1,ns
    read(22,*)nrcver(is)
    read(22,*)tmpstn,source(is)%x,source(is)%y,source(is)%z
    if((source(is)%x .gt. maxx) .or. &
      &(source(is)%x .lt. minx) .or. &
      &(source(is)%y .gt. maxy) .or. &
      &(source(is)%y .lt. miny) .or. &
      &(source(is)%z .gt. maxz) .or. &
      &(source(is)%z .lt. minz))then
        write(*,*)"Number ",is," source is outside study area!!"
        write(*,*)"Coordinates: ",source(is)
        stop
    end if
    do ir=1,nrcver(is)
        ndata=ndata+1
        read(22,*)tmpt,tmpx,tmpy,tmpz,evnid(ndata),receiver(ndata)%x,&
        &receiver(ndata)%y,receiver(ndata)%z
        evnum(ndata)=ndata
        if(evnid(ndata) .gt. ttlevn)then
            ttlevn=evnid(ndata)
        end if
        if((receiver(ndata)%x .gt. maxx) .or. &
          &(receiver(ndata)%x .lt. minx) .or. &
          &(receiver(ndata)%y .gt. maxy) .or. &
          &(receiver(ndata)%y .lt. miny) .or. &
          &(receiver(ndata)%z .gt. maxz) .or. &
          &(receiver(ndata)%z .lt. minz))then
            write(*,*)"Number ",ndata," receiver is outside study area!!"
            write(*,*)"Coordinates: ",receiver(ndata)
            stop
        end if
    end do
end do
close(22)

