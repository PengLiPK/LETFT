#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=syn_work
inpfdir=input_file
#######################################################

# Parameters of topography file
#######################################################
topo_file=hawaii_tp4tomo2.xyz
topo_minlon=203.80
topo_maxlon=205.383
topo_minlan=18.80
topo_maxlan=20.383
topo_xnum=96
topo_ynum=96
topo_depth=9.673
topo_vair=2.0
#######################################################

# Parameters of model evaluation
#######################################################
damp=100
nevn=2493
nsl=3264
ndata=25465
sta_num=35
minlon=204.50
maxlon=205.20
minlan=19.20
maxlan=19.70
minz=0
maxz=25.0
fmm_dx=0.005
fmm_dy=0.005
fmm_dz=0.5
fmm_rf_nx=5.0
fmm_rf_ny=5.0
fmm_rf_nz=5.0
fmm_rf_x=0.045
fmm_rf_y=0.045
fmm_rf_z=5
#######################################################


# Check working dir is exist or not.
if [ ! -d $wkdir ]
then
	echo "There is no dir called $wkdir!"
	echo "You need to run cp_runsyn_vel_t.bash first!"
	exit
fi
if [ ! -d $wkdir/syn_tomo/ ]
then
	echo "There is no dir called $wkdir/syn_tomo/!"
	echo "You need to run cp_runsyn_tomo.bash first!"
	exit
fi
		
# Copy files to syn_tomo dir.
cp ./src/mdresolmat ./$wkdir/syn_tomo/


# Enter work dir
cd ./$wkdir/syn_tomo/
pwd

# Model resolution matrix
nsldata=$(wc -l dtsep.txt | awk '{print $1}')
echo "dtsep.txt" > mdresolmat.inp
echo "metafd.dat" >> mdresolmat.inp
echo "fdcsm.dat" >> mdresolmat.inp
echo "mdres.txt" >> mdresolmat.inp
echo "$nsldata $nsl" >> mdresolmat.inp
echo "$damp" >> mdresolmat.inp

./mdresolmat
echo "Model resolution maxtrix calculation finished!!"

# Residual
echo "syndataloc.txt" > fmm_fw_regul2tprl.inp
echo "vel.txt" >> fmm_fw_regul2tprl.inp
echo "$sta_num" >> fmm_fw_regul2tprl.inp
echo "$fmm_dx $fmm_dy $fmm_dz" >> fmm_fw_regul2tprl.inp
echo "$fmm_rf_nx $fmm_rf_ny $fmm_rf_nz" >> fmm_fw_regul2tprl.inp
echo "$fmm_rf_x $fmm_rf_y $fmm_rf_z" >> fmm_fw_regul2tprl.inp
echo "$minlon $maxlon" >> fmm_fw_regul2tprl.inp
echo "$minlan $maxlan" >> fmm_fw_regul2tprl.inp
echo "$minz $maxz" >> fmm_fw_regul2tprl.inp
echo "2" >> fmm_fw_regul2tprl.inp
echo "2" >> fmm_fw_regul2tprl.inp
echo "$topo_file" >> fmm_fw_regul2tprl.inp
echo "$topo_minlon $topo_maxlon" >> fmm_fw_regul2tprl.inp
echo "$topo_minlan $topo_maxlan" >> fmm_fw_regul2tprl.inp
echo "$topo_xnum $topo_ynum" >> fmm_fw_regul2tprl.inp
echo "$topo_depth" >> fmm_fw_regul2tprl.inp
echo "$topo_vair" >> fmm_fw_regul2tprl.inp
./fmm_fw_regul2tprl
echo "Forward finished!"

awk '{if(NF==9)print $1,$9}' syndataloc.txt > tempt
paste t.txt tempt | awk '{print ($1+$6)-$5}' > res_final.txt

cp res_final.txt mdres.txt ../syn_outf/

# Exit work dir
cd ../../
