fc=gfortran

flags= -Wall -Wextra -Wconversion -c

objects1 = fmm_synt3d_v.o ../modules/strct.o ../modules/fmarch.o


all: fmm_synt3d_v


fmm_synt3d_v : $(objects1)
	$(fc) -o fmm_synt3d_v $(objects1)
fmm_synt3d_v.o : fmm_synt3d_v.f90
	$(fc) -c fmm_synt3d_v.f90
