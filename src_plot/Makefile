fc=gfortran

flags= -Wall -Wextra -Wconversion -c

objects1 = cutslice_reg2.o strct.o
objects2 = cutslice_reg2b.o strct.o
obj_crr2 = cut_ressl_reg2.o strct.o
objects4 = cutv_reg2.o strct.o
obj_cr1 = cutvdiff_reg1.o strct.o
obj_cr2 = cutvdiff_reg2.o strct.o
obj_cr1xyz = cutvdiff_reg1_xyz.o strct.o
obj_cr2xyz = cutvdiff_reg2_xyz.o strct.o
obj_crr2xyz = cut_ressl_reg2_xyz.o strct.o
obj_cl2xyz = cutvdiff_line2_xyz.o strct.o


all: strct.o cut_ressl_reg2 cutvdiff_reg1 cutvdiff_reg2 \
	cutvdiff_reg1_xyz cutvdiff_reg2_xyz cut_ressl_reg2_xyz \
	valmask cutvdiff_line2_xyz


cut_ressl_reg2 : $(obj_crr2)
	$(fc) -o cut_ressl_reg2 $(obj_crr2)
cut_ressl_reg2.o : cut_ressl_reg2.f90
	$(fc) $(flags) cut_ressl_reg2.f90
cutvdiff_reg1 : $(obj_cr1)
	$(fc) -o cutvdiff_reg1 $(obj_cr1)
cutvdiff_reg1.o : cutvdiff_reg1.f90
	$(fc) $(flags) cutvdiff_reg1.f90
cutvdiff_reg2 : $(obj_cr2)
	$(fc) -o cutvdiff_reg2 $(obj_cr2)
cutvdiff_reg2.o : cutvdiff_reg2.f90
	$(fc) $(flags) cutvdiff_reg2.f90
cutvdiff_reg1_xyz : $(obj_cr1xyz)
	$(fc) -o cutvdiff_reg1_xyz $(obj_cr1xyz)
cutvdiff_reg1_xyz.o : cutvdiff_reg1_xyz.f90
	$(fc) $(flags) cutvdiff_reg1_xyz.f90
cutvdiff_reg2_xyz : $(obj_cr2xyz)
	$(fc) -o cutvdiff_reg2_xyz $(obj_cr2xyz)
cutvdiff_reg2_xyz.o : cutvdiff_reg2_xyz.f90
	$(fc) $(flags) cutvdiff_reg2_xyz.f90
cut_ressl_reg2_xyz : $(obj_crr2xyz)
	$(fc) -o cut_ressl_reg2_xyz $(obj_crr2xyz)
cut_ressl_reg2_xyz.o : cut_ressl_reg2_xyz.f90
	$(fc) $(flags) cut_ressl_reg2_xyz.f90
valmask : valmask.o
	$(fc) -o valmask valmask.o
valmask.o : valmask.f90
	$(fc) $(flags) valmask.f90
cutvdiff_line2_xyz : $(obj_cl2xyz)
	$(fc) -o cutvdiff_line2_xyz $(obj_cl2xyz)
cutvdiff_line2_xyz.o : cutvdiff_line2_xyz.f90
	$(fc) $(flags) cutvdiff_line2_xyz.f90
strct.o : strct.f90
	$(fc) -c strct.f90
