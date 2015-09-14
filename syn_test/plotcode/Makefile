fc=gfortran

flags= -Wall -Wextra -Wconversion -c

objects1 = cutslice_reg2.o strct.o
objects2 = cutslice_reg2b.o strct.o
objects3 = cut_ressl_reg2.o strct.o
objects4 = cutv_reg2.o strct.o
objects5 = cutvdiff_reg1.o strct.o
objects6 = cutvdiff_reg2.o strct.o
objects7 = cutvdiff_reg1_xyz.o strct.o
objects8 = cutvdiff_reg2_xyz.o strct.o
objects9 = cut_ressl_reg2_xyz.o strct.o
objects10 = cutvdiff_line2_xyz.o strct.o


all: cutslice_reg2 cutslice_reg2b cut_ressl_reg2 cutv_reg2 cutvdiff_reg1 cutvdiff_reg2 cutvdiff_reg1_xyz cutvdiff_reg2_xyz cut_ressl_reg2_xyz valmask cutvdiff_line2_xyz


cutslice_reg2 : $(objects1)
	$(fc) -o cutslice_reg2 $(objects1)
cutslice_reg2.o : cutslice_reg2.f90
	$(fc) $(flags) cutslice_reg2.f90
cutslice_reg2b : $(objects2)
	$(fc) -o cutslice_reg2b $(objects2)
cutslice_reg2b.o : cutslice_reg2b.f90
	$(fc) $(flags) cutslice_reg2b.f90
cut_ressl_reg2 : $(objects3)
	$(fc) -o cut_ressl_reg2 $(objects3)
cut_ressl_reg2.o : cut_ressl_reg2.f90
	$(fc) $(flags) cut_ressl_reg2.f90
cutv_reg2 : $(objects4)
	$(fc) -o cutv_reg2 $(objects4)
cutv_reg2.o : cutv_reg2.f90
	$(fc) $(flags) cutv_reg2.f90
cutvdiff_reg1 : $(objects5)
	$(fc) -o cutvdiff_reg1 $(objects5)
cutvdiff_reg1.o : cutvdiff_reg1.f90
	$(fc) $(flags) cutvdiff_reg1.f90
cutvdiff_reg2 : $(objects6)
	$(fc) -o cutvdiff_reg2 $(objects6)
cutvdiff_reg2.o : cutvdiff_reg2.f90
	$(fc) $(flags) cutvdiff_reg2.f90
cutvdiff_reg1_xyz : $(objects7)
	$(fc) -o cutvdiff_reg1_xyz $(objects7)
cutvdiff_reg1_xyz.o : cutvdiff_reg1_xyz.f90
	$(fc) $(flags) cutvdiff_reg1_xyz.f90
cutvdiff_reg2_xyz : $(objects8)
	$(fc) -o cutvdiff_reg2_xyz $(objects8)
cutvdiff_reg2_xyz.o : cutvdiff_reg2_xyz.f90
	$(fc) $(flags) cutvdiff_reg2_xyz.f90
cut_ressl_reg2_xyz : $(objects9)
	$(fc) -o cut_ressl_reg2_xyz $(objects9)
cut_ressl_reg2_xyz.o : cut_ressl_reg2_xyz.f90
	$(fc) $(flags) cut_ressl_reg2_xyz.f90
valmask : valmask.o
	$(fc) -o valmask valmask.o
valmask.o : valmask.f90
	$(fc) $(flags) valmask.f90
cutvdiff_line2_xyz : $(objects10)
	$(fc) -o cutvdiff_line2_xyz $(objects10)
cutvdiff_line2_xyz.o : cutvdiff_line2_xyz.f90
	$(fc) $(flags) cutvdiff_line2_xyz.f90
strct.o : strct.f90
	$(fc) -c strct.f90
