program testtriint

use strct
implicit none
real(kind=8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
real(kind=8) :: x,y,z



open(21,file='testtriint.inp',status='old')
read(21,*)x1,y1,z1
read(21,*)x2,y2,z2
read(21,*)x3,y3,z3
read(21,*)x,y
close(21)


z=trisurfintp(z1,z2,z3,x1,y1,x2,y2,x3,y3,x,y)

print *, z


stop
end
