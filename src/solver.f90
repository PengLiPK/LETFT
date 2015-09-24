module solver

integer,parameter :: maxrow=60
integer,parameter :: maxcol=4

contains

! QR decomposition by Householder transform.
! Return qt, r is stored in a.
! a: m by n and m > n
!--------------------------------------------------------------------------
subroutine qr(qt,a,nrow,ncol)

implicit none
real(kind=8),intent(inout) :: qt(maxrow,maxrow),a(maxrow,maxcol)
real(kind=8),allocatable :: c(:,:),h(:,:)
real(kind=8),allocatable :: tmpvec(:)
real(kind=8) :: xp
integer,intent(in) :: nrow,ncol
integer :: i,j,k

allocate(c(nrow,nrow))
allocate(h(nrow,nrow))
allocate(tmpvec(nrow))

tmpvec(1:nrow)=a(1:nrow,1)
call householder(xp,h,tmpvec,nrow)
do i=1,nrow
    do j=1,nrow
        qt(i,j)=h(i,j)
    end do
end do
!qt(1:nrow,1:nrow)=h(1:nrow,1:nrow)
a(1,1)=xp
a(2:nrow,1)=0d0
c=matmul(qt(1:nrow,1:nrow),a(1:nrow,2:ncol))
do i=1,nrow
    do j=2,ncol
        a(i,j)=c(i,j-1)
    end do
end do
!a(1:nrow,2:ncol)=c(1:nrow,1:ncol-1)

do k=2,ncol
    do i=1,nrow+1-k
        tmpvec(i)=a(i+k-1,k)
    end do
    !tmpvec(1:nrow+1-k)=a(k:nrow,k)
    call householder(xp,h,tmpvec,nrow+1-k)
    
    c=matmul(h(1:nrow+1-k,1:nrow+1-k),qt(k:nrow,1:nrow))
    do i=k,nrow
        do j=1,nrow
            qt(i,j)=c(i-k+1,j)
        end do
    end do
    !qt(k:nrow,1:nrow)=c(1:nrow+1-k,1:nrow)

    a(k,k)=xp
    a(k+1:nrow,k)=0d0
    if(k .lt. ncol)then
        c=matmul(h(1:nrow+1-k,1:nrow+1-k),a(k:nrow,k+1:ncol))
        do i=k,nrow
            do j=k+1,ncol
                a(i,j)=c(i-k+1,j-k)
            end do
        end do
        !a(k:nrow,k+1:ncol)=c(1:nrow+1-k,1:ncol-k)
    end if
end do

deallocate(c)
deallocate(h)
deallocate(tmpvec)
return
end subroutine
!--------------------------------------------------------------------------


! Euclidean Norm of a vector
!--------------------------------------------------------------------------
subroutine householder(xp,h,x,len)

implicit none
real(kind=8),allocatable,intent(in) :: x(:)
real(kind=8),allocatable,intent(inout) :: h(:,:)
real(kind=8),allocatable :: u(:)
real(kind=8) :: xp
integer,intent(in) :: len
integer :: i,j


allocate(u(len))

! Calculate u.
if(x(1) .ge. 0d0)then
    xp=-vec_norm2(x,len)
else
    xp=vec_norm2(x,len)
end if

u(1)=sqrt((xp-x(1))/(2.0d0*xp))

do i=2,len
    u(i)=-x(i)/(2.0d0*xp*u(1))
end do

! H=I-2uu'
do i=1,len
    do j=1,len
        if(i .eq. j)then
            h(i,j)=1.0d0-2.0d0*u(i)*u(j)
        else
            h(i,j)=-2.0d0*u(i)*u(j)
        end if
    end do
end do


deallocate(u)

return
end subroutine householder
!--------------------------------------------------------------------------



! Matrix multiply, results are stored in c.
!--------------------------------------------------------------------------
subroutine matrix_mul(c,a,arow,acol,b,brow,bcol)

implicit none
real(kind=8),allocatable,intent(in):: a(:,:),b(:,:)
real(kind=8),allocatable,intent(inout) :: c(:,:)
integer,intent(in) :: arow,acol,brow,bcol
integer :: i,j,k

if(acol .ne. brow)then
    print *, "Error. Col of matrix a != Row of matrix b"
    print *, acol,brow
    stop
end if

do i=1,arow
    do j=1,bcol
        c(i,j)=0d0
        do k=1,acol
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
    end do
end do

return
end subroutine matrix_mul
!--------------------------------------------------------------------------



! Euclidean Norm of a vector
!--------------------------------------------------------------------------
real(kind=8) function vec_norm2(x,len)

implicit none
real(kind=8),allocatable,intent(in) :: x(:)
integer,intent(in) :: len
integer :: i


vec_norm2=0d0
do i=1,len
    vec_norm2=vec_norm2+x(i)*x(i)
end do

vec_norm2=sqrt(vec_norm2)

return
end function vec_norm2
!--------------------------------------------------------------------------


! Inner product of two vectors
!--------------------------------------------------------------------------
real(kind=8) function vec_in_prod(x,y,len)

implicit none
real(kind=8),allocatable,intent(in) :: x(:),y(:)
integer,intent(in) :: len
integer :: i


vec_in_prod=0d0
do i=1,len
    vec_in_prod=vec_in_prod+x(i)*y(i)
end do

return
end function vec_in_prod
!--------------------------------------------------------------------------


end module solver
!--------------------------------------------------------------------------
