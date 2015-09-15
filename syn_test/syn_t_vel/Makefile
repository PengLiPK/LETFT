fc=gfortran

flags= -Wall -Wextra -Wconversion -c

objects1 = fmm_synt3d_v.o strct.o fmarch.o
objects2 = prevel3d.o strct.o


all: fmm_synt3d_v prevel3d


fmm_synt3d_v : $(objects1)
	$(fc) -o fmm_synt3d_v $(objects1)
fmm_synt3d_v.o : fmm_synt3d_v.f90
	$(fc) -c fmm_synt3d_v.f90
prevel3d : $(objects2)
	$(fc) -o prevel3d $(objects2)
prevel3d.o : prevel3d.f90
	$(fc) -c prevel3d.f90

