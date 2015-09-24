#!/bin/bash

time_begin=$(date +%s.%N)

# Parameters of working dir
#######################################################
wkdir=syn_work_ps
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
topo_vair_p=2.0
topo_vair_s=2.0
#######################################################

# Parameters of model evaluation
#######################################################
dampp=200
damps=150
nevn=1800
npmodel=3264
sta_num_p=35
sta_num_s=34
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
cp ./src/mdresolmat_ps ./$wkdir/syn_tomo/


# Enter work dir
cd ./$wkdir/syn_tomo/
pwd

# Model resolution matrix
nsldata=$(wc -l dtsep.txt | awk '{print $1}')
echo "dtsep.txt" > mdresolmat_ps.inp
echo "metafd.dat" >> mdresolmat_ps.inp
echo "fdcsm.dat" >> mdresolmat_ps.inp
echo "mdres.txt" >> mdresolmat_ps.inp
echo "$nsldata $npmodel $npmodel" >> mdresolmat_ps.inp
echo "$dampp $damps" >> mdresolmat_ps.inp

./mdresolmat_ps
head -$npmodel mdres.txt > mdres_p.txt
tail -$npmodel mdres.txt > mdres_s.txt
echo "Model resolution maxtrix calculation finished!!"

# Residual
echo "syndataloc_p.txt" > fmm_fw_regul2tprl.inp
echo "vel_p.txt" >> fmm_fw_regul2tprl.inp
echo "$sta_num_p" >> fmm_fw_regul2tprl.inp
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
echo "$topo_vair_p" >> fmm_fw_regul2tprl.inp
./fmm_fw_regul2tprl
mv fd.dat fd_p.dat
mv t.txt t_p.txt

echo "syndataloc_s.txt" > fmm_fw_regul2tprl.inp
echo "vel_s.txt" >> fmm_fw_regul2tprl.inp
echo "$sta_num_s" >> fmm_fw_regul2tprl.inp
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
echo "$topo_vair_s" >> fmm_fw_regul2tprl.inp
./fmm_fw_regul2tprl
mv fd.dat fd_s.dat
mv t.txt t_s.txt

echo "Forward finished!"


cat syndataloc_p.txt syndataloc_s.txt |  awk '{if(NF==9)print $1,$9}' > tempt
cat t_p.txt t_s.txt | awk '{print $1}' > t.txt
paste t.txt tempt | awk '{print ($1+$3)-$2}' > res_final.txt

cp res_final.txt mdres_p.txt mdres_s.txt ../syn_outf/

# Exit work dir
cd ../../

time_end=$(date +%s.%N)
time_used=$(echo "scale=5;($time_end - $time_begin)/60.0" | bc)
echo "Total time used: $time_used min."
