fc=gfortran
libdir=/nethome/pli/lib

obj_fmmsynt3d = fmm_synt3d_v.o strct.o fmarch.o
obj_prevel3d = prevel3d.o strct.o
obj_genGSnoise = genGSnoise.o strct.o
obj_migraloc = migraloc.o strct.o
obj_fmmfwregul2tprl = fmm_fw_regul2tprl.o strct.o fmarch.o
obj_invctrlds = invctrlds.o strct.o
obj_dloc2locmulti = dloc2loc_multi.o strct.o


all: strct.o fmarch.o solver.o fmm_synt3d_v prevel3d genGSnoise\
	 migraloc fmm_fw_regul2tprl fd2csm_row fdsep_qr fdsep_qr_wei\
	 fd2csm calmdres2 invsvddt_damp invsvddt_damp_ps \
	 invctrlds invloc_svd_wei dloc2loc_multi mdresolmat \
	 mdresolmat_ps comb_ps initinpv

fmm_synt3d_v : $(obj_fmmsynt3d)
	$(fc) -o fmm_synt3d_v $(obj_fmmsynt3d)
fmm_synt3d_v.o : fmm_synt3d_v.f90
	$(fc) -c fmm_synt3d_v.f90
prevel3d : $(obj_prevel3d)
	$(fc) -o prevel3d $(obj_prevel3d)
prevel3d.o : prevel3d.f90
	$(fc) -c prevel3d.f90
genGSnoise : $(obj_genGSnoise)
	$(fc) -o genGSnoise $(obj_genGSnoise)
genGSnoise.o : genGSnoise.f90
	$(fc) -c genGSnoise.f90
migraloc : $(obj_migraloc)
	$(fc) -o migraloc $(obj_migraloc)
migraloc.o : migraloc.f90
	$(fc) -c migraloc.f90
fmm_fw_regul2tprl : $(obj_fmmfwregul2tprl)
	$(fc) -o fmm_fw_regul2tprl $(obj_fmmfwregul2tprl)
fmm_fw_regul2tprl.o : fmm_fw_regul2tprl.f90
	$(fc) -c fmm_fw_regul2tprl.f90
fd2csm_row : fd2csm_row.o
	$(fc) -o fd2csm_row fd2csm_row.o
fd2csm_row.o : fd2csm_row.f90
	$(fc) -c fd2csm_row.f90
fdsep_qr : fdsep_qr.f90 solver.o
	$(fc) -mcmodel=medium -o fdsep_qr fdsep_qr.f90 solver.o
fdsep_qr_wei : fdsep_qr_wei.f90 solver.o strct.o
	$(fc) -o fdsep_qr_wei fdsep_qr_wei.f90 solver.o strct.o
fd2csm : fd2csm.o
	$(fc) -o fd2csm fd2csm.o
fd2csm.o : fd2csm.f90
	$(fc) -c fd2csm.f90
calmdres2 : calmdres2.o
	$(fc) -o calmdres2 calmdres2.o
calmdres2.o : calmdres2.f90
	$(fc) -c calmdres2.f90
invsvddt_damp : invsvddt_damp.f90
	$(fc) -mcmodel=medium -L$(libdir) -o invsvddt_damp invsvddt_damp.f90 -llapack -lrefblas
invsvddt_damp_ps : invsvddt_damp_ps.f90
	$(fc) -mcmodel=medium -L$(libdir) -o invsvddt_damp_ps invsvddt_damp_ps.f90 -llapack -lrefblas
invctrlds : $(obj_invctrlds)
	$(fc) -o invctrlds $(obj_invctrlds)
invctrlds.o : invctrlds.f90
	$(fc) -c invctrlds.f90
invloc_svd_wei : invloc_svd_wei.f90 strct.o
	$(fc) -mcmodel=medium -L$(libdir) -o invloc_svd_wei invloc_svd_wei.f90 strct.o -llapack -lrefblas
dloc2loc_multi : $(obj_dloc2locmulti)
	$(fc) -o dloc2loc_multi $(obj_dloc2locmulti)
dloc2loc_multi.o : dloc2loc_multi.f90
	$(fc) -c dloc2loc_multi.f90
mdresolmat : mdresolmat.f90
	$(fc) -mcmodel=medium -L$(libdir) -o mdresolmat mdresolmat.f90 -llapack -lrefblas
mdresolmat_ps : mdresolmat_ps.f90
	$(fc) -mcmodel=medium -L$(libdir) -o mdresolmat_ps mdresolmat_ps.f90 -llapack -lrefblas
comb_ps : comb_ps.f90
	$(fc) -o comb_ps comb_ps.f90
initinpv : initinpv.f90 strct.o
	$(fc) -o initinpv initinpv.f90 strct.o
strct.o : strct.f90
	$(fc) -c strct.f90
fmarch.o : fmarch.f90
	$(fc) -c fmarch.f90
solver.o : solver.f90
	$(fc) -c solver.f90
