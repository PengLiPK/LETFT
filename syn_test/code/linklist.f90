module linklist
    implicit none
    type :: node_linklist
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: vel
        type(node_linklist),pointer :: prev
        type(node_linklist),pointer :: next
    end type node_linklist

    ! 3D nodes.
    type :: node_linklist3d
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: z
        real(kind=8) :: vel
        type(node_linklist3d),pointer :: prev
        type(node_linklist3d),pointer :: next
    end type node_linklist3d
    
contains

subroutine list_output(list,fname)
    implicit none
    type(node_linklist),pointer :: list,p
    character(len=70) :: fname
    
    if(associated(list%prev) .eqv. .false.)then
        p=>list%next
    else
        p=>list
    end if
    
    if(trim(fname) .eq. 'screen')then
        do while(associated(p%next))
            write(*,*)p%x,p%y,p%vel
            p=>p%next
        end do
    else
        open(900,file=fname)
        do while(associated(p%next))
            write(900,*)p%x,p%y,p%vel
            p=>p%next
        end do
    end if
    
    return

end subroutine list_output


! Delete
subroutine list_del(item)
    implicit none
    type(node_linklist),pointer :: item
    type(node_linklist),pointer :: prev,next
    
    prev=>item%prev
    next=>item%next
    deallocate(item)
    if(associated(prev))prev%next=>next
    if(associated(next))next%prev=>prev
    item=>next

    return

end subroutine list_del

! Insert item behind p
subroutine list_ins(p,x,y,vel,bfswitch)
    implicit none
    type(node_linklist),pointer :: p
    type(node_linklist),pointer :: newnode
    real(kind=8) :: x,y,vel
    integer :: errormem
    integer :: bfswitch


    allocate(newnode,stat=errormem)
    if(errormem/=0)then
        write(*,*)"Out of memory!! (linklist.f90)"
        stop
    end if
    newnode%x=x
    newnode%y=y
    newnode%vel=vel
    
    if(bfswitch .eq. 1)then
        nullify(newnode%next)
        if(associated(p%next))then
            newnode%next=>p%next
            p%next%prev=>newnode
        end if
        newnode%prev=>p
        p%next=>newnode
    else if(bfswitch .eq. 2)then
        nullify(newnode%prev)
        if(associated(p%prev))then
            newnode%prev=>p%prev
            p%prev%next=>newnode
        end if
        newnode%next=>p
        p%prev=>newnode
    end if


    return

end subroutine list_ins

! Find max element by x in list
function list_max(listhead)
    implicit none
    type(node_linklist),pointer :: listhead
    type(node_linklist),pointer :: list_max

    list_max=>listhead
    do while(associated(list_max%next))
        if(list_max%next%x .gt. list_max%x)then
            list_max=>list_max%next
        end if
    end do

    return

end function list_max



! For 3D nodes
subroutine list_output3d(list,fname)
    implicit none
    type(node_linklist3d),pointer :: list,p
    character(len=70) :: fname
    
    if(associated(list%prev) .eqv. .false.)then
        p=>list%next
    else
        p=>list
    end if
    
    if(trim(fname) .eq. 'screen')then
        do while(associated(p%next))
            write(*,*)p%x,p%y,p%z,p%vel
            p=>p%next
        end do
    else
        open(900,file=fname)
        do while(associated(p%next))
            write(900,*)p%x,p%y,p%z,p%vel
            p=>p%next
        end do
    end if
    
    return

end subroutine list_output3d


! Delete
subroutine list_del3d(item)
    implicit none
    type(node_linklist3d),pointer :: item
    type(node_linklist3d),pointer :: prev,next
    
    prev=>item%prev
    next=>item%next
    deallocate(item)
    if(associated(prev))prev%next=>next
    if(associated(next))next%prev=>prev
    item=>next

    return

end subroutine list_del3d

! Insert item behind p
subroutine list_ins3d(p,x,y,z,vel,bfswitch)
    implicit none
    type(node_linklist3d),pointer :: p
    type(node_linklist3d),pointer :: newnode
    real(kind=8) :: x,y,z,vel
    integer :: errormem
    integer :: bfswitch


    allocate(newnode,stat=errormem)
    if(errormem/=0)then
        write(*,*)"Out of memory!! (linklist.f90)"
        stop
    end if
    newnode%x=x
    newnode%y=y
    newnode%z=z
    newnode%vel=vel
    
    if(bfswitch .eq. 1)then
        nullify(newnode%next)
        if(associated(p%next))then
            newnode%next=>p%next
            p%next%prev=>newnode
        end if
        newnode%prev=>p
        p%next=>newnode
    else if(bfswitch .eq. 2)then
        nullify(newnode%prev)
        if(associated(p%prev))then
            newnode%prev=>p%prev
            p%prev%next=>newnode
        end if
        newnode%next=>p
        p%prev=>newnode
    end if


    return

end subroutine list_ins3d

! Find max element by x in list
function list_max3d(listhead)
    implicit none
    type(node_linklist3d),pointer :: listhead
    type(node_linklist3d),pointer :: list_max3d

    list_max3d=>listhead
    do while(associated(list_max3d%next))
        if(list_max3d%next%x .gt. list_max3d%x)then
            list_max3d=>list_max3d%next
        end if
    end do

    return

end function list_max3d


end module



