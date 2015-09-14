fc=gfortran
#flags= -Wall -Wextra -Wconversion -c

objects2 = regul_grid.o strct.o
objects3 = fmm_fw_regul.o strct.o fmarch.o
objects4 = del0vol.o strct.o
objects5 = fmm_fw_regul2.o strct.o fmarch.o
objects6 = fmm_fw_regul3.o strct.o fmarch.o
objects7 = fmm_fw_regul2tp.o strct.o fmarch.o
objects8 = fmm_res.o strct.o fmarch.o
objects9 = invctrl.o strct.o
objects10 = initinpv.o strct.o
objects11 = fmm_fw_regul2tprl.o strct.o fmarch.o
objects12 = dloc2loc.o strct.o
objects13 = dloc2loc_multi.o strct.o
objects14 = invctrlds.o strct.o
objects15 = genGSnoise.o strct.o
objects16 = migraloc.o strct.o
objects17 = fmm_fw_regul2rl.o strct.o fmarch.o
objects18 = dloc2loc_multi_notp.o strct.o


#all: fmm_adpt3d fmm_fw_regul del0vol inverseprb
all: regul_grid fmm_fw_regul del0vol inverseprb inverseprbdt invpcgdt fd2csm cal_fdx invqrdt invsvddt fmm_fw_regul2 fmm_fw_regul2tp dv2ds fmm_fw_regul3 cal_fds_lsres cal_fds ds2whvel calmdres fmm_res readdata invctrl calmdres2 initinpv cutpro4inv cvtinvout fmm_fw_regul2tprl fd2locfd dloc2loc dloc2loc_multi invctrlds genGSnoise addweight invcgls migraloc invcgls2 fd2csm_row fdsep_svd fdsep_qr fdsep_qr_wei fdsplit invloc invcgls3 invloc_qr invloc_qr_wei invloc_svd invloc_svd_wei invsvddt_damp invsvddt2 invsvddt2_damp mdresolmat fmm_fw_regul2rl dloc2loc_multi_notp addns

regul_grid: $(objects2)
	$(fc) -o regul_grid $(objects2)
regul_grid.o: regul_grid.f90
	$(fc) -c regul_grid.f90
fmm_fw_regul : $(objects3)
	$(fc) -o fmm_fw_regul $(objects3)
fmm_fw_regul.o : fmm_fw_regul.f90
	$(fc) -c fmm_fw_regul.f90
strct.o : strct.f90
	$(fc) -c strct.f90
linklist.o : linklist.f90
	$(fc) -c linklist.f90
fmarch.o : fmarch.f90
	$(fc) -c fmarch.f90
del0vol : $(objects4)
	$(fc) -o del0vol $(objects4)
del0vol.o : del0vol.f90
	$(fc) -c del0vol.f90
inverseprb : inverseprb.o
	$(fc) -o inverseprb inverseprb.o
inverseprb.o : inverseprb.f90
	$(fc) -c inverseprb.f90
inverseprbdt : inverseprbdt.o
	$(fc) -o inverseprbdt inverseprbdt.o
inverseprbdt.o : inverseprbdt.f90
	$(fc) -c inverseprbdt.f90
invpcgdt : invpcgdt.o
	$(fc) -o invpcgdt invpcgdt.o
invpcgdt.o : invpcgdt.f90
	$(fc) -c invpcgdt.f90
fd2csm : fd2csm.o
	$(fc) -o fd2csm fd2csm.o
fd2csm.o : fd2csm.f90
	$(fc) -c fd2csm.f90
cal_fdx : cal_fdx.o
	$(fc) -o cal_fdx cal_fdx.o
cal_fdx.o : cal_fdx.f90
	$(fc) -c cal_fdx.f90
invqrdt : invqrdt.o
	$(fc) -L/nethome/pli/lib -o invqrdt invqrdt.o -llapack -lrefblas
invqrdt.o : invqrdt.f90
	$(fc) -c invqrdt.f90
invsvddt : invsvddt.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invsvddt invsvddt.f90 -llapack -lrefblas
fmm_fw_regul2 : $(objects5)
	$(fc) -o fmm_fw_regul2 $(objects5)
fmm_fw_regul2.o : fmm_fw_regul2.f90
	$(fc) -c fmm_fw_regul2.f90
fmm_fw_regul2tp : $(objects7)
	$(fc) -o fmm_fw_regul2tp $(objects7)
fmm_fw_regul2tp.o : fmm_fw_regul2tp.f90
	$(fc) -c fmm_fw_regul2tp.f90
dv2ds : dv2ds.o
	$(fc) -o dv2ds dv2ds.o strct.o
dv2ds.o : dv2ds.f90
	$(fc) -c dv2ds.f90
fmm_fw_regul3 : $(objects6)
	$(fc) -o fmm_fw_regul3 $(objects6)
fmm_fw_regul3.o : fmm_fw_regul3.f90
	$(fc) -c fmm_fw_regul3.f90
cal_fds_lsres : cal_fds_lsres.o
	$(fc) -o cal_fds_lsres cal_fds_lsres.o
cal_fds_lsres.o : cal_fds_lsres.f90
	$(fc) -c cal_fds_lsres.f90
cal_fds : cal_fds.o
	$(fc) -o cal_fds cal_fds.o
cal_fds.o : cal_fds.f90
	$(fc) -c cal_fds.f90
ds2whvel : ds2whvel.o strct.o
	$(fc) -o ds2whvel ds2whvel.o strct.o
ds2whvel.o : ds2whvel.f90
	$(fc) -c ds2whvel.f90
calmdres : calmdres.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o calmdres calmdres.f90 -llapack -lrefblas
fmm_res : $(objects8)
	$(fc) -o fmm_res $(objects8)
fmm_res.o : fmm_res.f90
	$(fc) -c fmm_res.f90
readdata : readdata.o
	$(fc) -o readdata readdata.o
readdata.o : readdata.f90
	$(fc) -c readdata.f90
invctrl : $(objects9)
	$(fc) -o invctrl $(objects9)
invctrl.o : invctrl.f90
	$(fc) -c invctrl.f90
calmdres2 : calmdres2.o
	$(fc) -o calmdres2 calmdres2.o
calmdres2.o : calmdres2.f90
	$(fc) -c calmdres2.f90
initinpv : $(objects10)
	$(fc) -o initinpv $(objects10)
initinpv.o : initinpv.f90
	$(fc) -c initinpv.f90
cutpro4inv : cutpro4inv.o
	$(fc) -o cutpro4inv cutpro4inv.o
cutpro4inv.o : cutpro4inv.f90
	$(fc) -c cutpro4inv.f90
cvtinvout : cvtinvout.o
	$(fc) -o cvtinvout cvtinvout.o
cvtinvout.o : cvtinvout.f90
	$(fc) -c cvtinvout.f90
fmm_fw_regul2tprl : $(objects11)
	$(fc) -o fmm_fw_regul2tprl $(objects11)
fmm_fw_regul2tprl.o : fmm_fw_regul2tprl.f90
	$(fc) -c fmm_fw_regul2tprl.f90
fd2locfd : fd2locfd.o
	$(fc) -o fd2locfd fd2locfd.o
fd2locfd.o : fd2locfd.f90
	$(fc) -c fd2locfd.f90
dloc2loc : $(objects12)
	$(fc) -o dloc2loc $(objects12)
dloc2loc.o : dloc2loc.f90
	$(fc) -c dloc2loc.f90
dloc2loc_multi : $(objects13)
	$(fc) -o dloc2loc_multi $(objects13)
dloc2loc_multi.o : dloc2loc_multi.f90
	$(fc) -c dloc2loc_multi.f90
invctrlds : $(objects14)
	$(fc) -o invctrlds $(objects14)
invctrlds.o : invctrlds.f90
	$(fc) -c invctrlds.f90
genGSnoise : $(objects15)
	$(fc) -o genGSnoise $(objects15)
genGSnoise.o : genGSnoise.f90
	$(fc) -c genGSnoise.f90
addweight : addweight.o
	$(fc) -o addweight addweight.o
addweight.o : addweight.f90
	$(fc) -c addweight.f90
migraloc : $(objects16)
		$(fc) -o migraloc $(objects16)
migraloc.o : migraloc.f90
		$(fc) -c migraloc.f90
invcgls : invcgls.o
		$(fc) -o invcgls invcgls.o
invcgls.o : invcgls.f90
		$(fc) -c invcgls.f90
invcgls2 : invcgls2.o
		$(fc) -o invcgls2 invcgls2.o
invcgls2.o : invcgls2.f90
		$(fc) -c invcgls2.f90
fd2csm_row : fd2csm_row.o
		$(fc) -o fd2csm_row fd2csm_row.o
fd2csm_row.o : fd2csm_row.f90
		$(fc) -c fd2csm_row.f90
fdsep_svd : fdsep_svd.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o fdsep_svd fdsep_svd.f90 -llapack -lrefblas
fdsep_qr : fdsep_qr.f90 solver.o
	$(fc) -mcmodel=medium -o fdsep_qr fdsep_qr.f90 solver.o
fdsep_qr_wei : fdsep_qr_wei.f90 solver.o strct.o
	$(fc) -mcmodel=medium -o fdsep_qr_wei fdsep_qr_wei.f90 solver.o strct.o
solver.o : solver.f90
		$(fc) -c solver.f90
fdsplit : fdsplit.o
		$(fc) -o fdsplit fdsplit.o
fdsplit.o : fdsplit.f90
		$(fc) -c fdsplit.f90
invloc : invloc.o
		$(fc) -o invloc invloc.o
invloc.o : invloc.f90
		$(fc) -c invloc.f90
invcgls3 : invcgls3.o
		$(fc) -o invcgls3 invcgls3.o
invcgls3.o : invcgls3.f90
		$(fc) -c invcgls3.f90
invloc_qr : invloc_qr.f90 solver.o
	$(fc) -mcmodel=medium -o invloc_qr invloc_qr.f90 solver.o
invloc_qr_wei : invloc_qr_wei.f90 solver.o strct.o
	$(fc) -mcmodel=medium -o invloc_qr_wei invloc_qr_wei.f90 solver.o strct.o
invloc_svd : invloc_svd.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invloc_svd invloc_svd.f90 -llapack -lrefblas
invloc_svd_wei : invloc_svd_wei.f90 strct.o
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invloc_svd_wei invloc_svd_wei.f90 strct.o -llapack -lrefblas
invsvddt_damp : invsvddt_damp.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invsvddt_damp invsvddt_damp.f90 -llapack -lrefblas
invsvddt2 : invsvddt2.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invsvddt2 invsvddt2.f90 -llapack -lrefblas
invsvddt2_damp : invsvddt2_damp.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o invsvddt2_damp invsvddt2_damp.f90 -llapack -lrefblas
mdresolmat : mdresolmat.f90
	$(fc) -mcmodel=medium -L/nethome/pli/lib -o mdresolmat mdresolmat.f90 -llapack -lrefblas
fmm_fw_regul2rl : $(objects17)
	$(fc) -o fmm_fw_regul2rl $(objects17)
fmm_fw_regul2rl.o : fmm_fw_regul2rl.f90
	$(fc) -c fmm_fw_regul2rl.f90
dloc2loc_multi_notp : $(objects18)
	$(fc) -o dloc2loc_multi_notp $(objects18)
dloc2loc_multi_notp.o : dloc2loc_multi_notp.f90
	$(fc) -c dloc2loc_multi_notp.f90
addns : addns.f90
		$(fc) -o addns addns.f90
