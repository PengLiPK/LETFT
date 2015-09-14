! Order fd.
program fdcsord


implicit none




! Input parameters
!----------------------------------------------------

open(22,file='invcgls.inp')
read(22,*)metafdfile
read(22,*)fdfile
read(22,*)vfdfile
read(22,*)nevn
close(22)


! Read fd
metag=0
open(40,file=metagf,status='old')
imt=0
izero=0
do im=1,nevn*4
    read(40,*)tmpmetag
    if(tmpmetag .ne. 0)then
        imt=imt+1
        metag(imt)=tmpmetag
    else
        izero=izero+1
        metagzero(izero)=im
    end if
end do
imtlocmax=imt
izerolocmax=izero
do im=nevn*4+1,immax
    read(40,*)tmpmetag
    if(tmpmetag .ne. 0)then
        imt=imt+1
        metag(imt)=tmpmetag
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
i=0
do imt=1,imtlocmax/4
    do icol=1,4
        do irow=1,metag((imt-1)*4+icol)
            i=i+1
            read(43,rec=i,iostat=status1)gtmp,ir
            if(status1/=0)exit
            gindex(ir,1)=irow
            gindex(ir,2)=imt
            g(irow,icol,imt)=gtmp
        end do
    end do
end do
icol=4
do imt=imtlocmax+1,imtmax
    icol=icol+1
    do irow=1,metag(imt)
        i=i+1
        read(43,rec=i,iostat=status1)gtmp,ir
        if(status1/=0)exit
        g(gindex(ir,1),icol,gindex(ir,2))=gtmp
    end do
end do

close(43)




do imt=imlocmax/4
    call dgesvd('S','A',metag((imt-1)*4+1),4,g(:,1:4,imt),lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)

    do im
