#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work_p
inpfdir=input_file
#######################################################

# Parameters of input file
######################################################
initvel_file=init_vel.txt
data_file=p_data.txt
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

# Parameters of tomography
#######################################################
damp=200
nevn=1800
nsl=3264
ndata=30102
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
	echo "You need to run cp_run_init.bash first!"
	exit
fi
if [ ! -d $wkdir/tomo/ ]
then
	mkdir $wkdir/tomo/
fi
		
# Copy input files
cp ./$inpfdir/$topo_file ./$wkdir/tomo/
cp ./$wkdir/outf/$initvel_file ./$wkdir/tomo/
cp ./$wkdir/outf/$data_file ./$wkdir/tomo/


# Copy files to tomography directory
cp ./src/fmm_fw_regul2tprl ./$wkdir/tomo/
cp ./src/fd2csm_row ./$wkdir/tomo/
cp ./src/fdsep_qr ./$wkdir/tomo/
cp ./src/fd2csm ./$wkdir/tomo/
cp ./src/calmdres2 ./$wkdir/tomo/
cp ./src/invsvddt_damp ./$wkdir/tomo/
cp ./src/invctrlds ./$wkdir/tomo/
cp ./src/invloc_svd_wei ./$wkdir/tomo/
cp ./src/dloc2loc_multi ./$wkdir/tomo/


# Run tomography program
# Enter work dir
cd ./$wkdir/tomo/
pwd
mv $data_file syndataloc.txt
mv $initvel_file vel.txt

cp syndataloc.txt syndataloc_0.txt
for((i=1;i<=2;i++))
do

# Forward
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

# Separate slowness derivative matrix form fd

	awk '{if(NF==9)print $1,$9}' syndataloc.txt > tempt
	paste t.txt tempt | awk '{print ($1+$6)-$5}' > dt.txt
	j=$(( $i - 1 ))
	cp dt.txt res_$j.txt

	echo "fd.dat" > fd2csm_row.inp
	echo "metafd_row.dat" >> fd2csm_row.inp
	echo "fdcsm_row.dat" >> fd2csm_row.inp
	echo "$ndata" >> fd2csm_row.inp
	./fd2csm_row

	
	echo "metafd_row.dat" > fdsep_qr.inp
	echo "fdcsm_row.dat" >> fdsep_qr.inp
	echo "dt.txt" >> fdsep_qr.inp
	echo "fdloc.dat" >> fdsep_qr.inp
	echo "dtloc.txt" >> fdsep_qr.inp
	echo "rloc.dat" >> fdsep_qr.inp
	echo "fdsep.dat" >> fdsep_qr.inp
	echo "dtsep.txt" >> fdsep_qr.inp
	echo "$ndata $nevn $nsl" >> fdsep_qr.inp
	./fdsep_qr
   

# Prepare input for inverse

	echo "fdsep.dat" > fd2csm.inp
	echo "metafd.dat" >> fd2csm.inp
	echo "fdcsm.dat" >> fd2csm.inp
	echo "$nsl" >> fd2csm.inp
	./fd2csm


	nsldata=$(wc -l dtsep.txt | awk '{print $1}')

# Run calmdres2, calculating DWS values, ray counts
    echo "metafd.dat" > calmdres2.inp
    echo "fdcsm.dat" >> calmdres2.inp
    echo "rct.txt" >> calmdres2.inp
    echo "dws.txt" >> calmdres2.inp
    echo "$nsldata $nsl" >> calmdres2.inp
  
    ./calmdres2
    echo "Run calmdres2 finished!!"

	echo " Prepare input for inverse finished!"


# Run invsvddt_damp
    echo "dtsep.txt" > invsvddt_damp.inp
    echo "metafd.dat" >> invsvddt_damp.inp
	echo "fdcsm.dat" >> invsvddt_damp.inp
    echo "out_dsl.txt" >> invsvddt_damp.inp
    echo "$nsldata $nsl" >> invsvddt_damp.inp
    echo "$damp" >> invsvddt_damp.inp

    ./invsvddt_damp


	echo "Inversion part finished!!"
	
	cp out_dsl.txt out_dsl$i.txt
	cp vel.txt vel_$i.txt

# Convert results of invsvddt to real velocity
    echo "vel.txt" > invctrlds.inp
	echo "out_dsl.txt" >> invctrlds.inp
	echo "vel_$i.txt" >> invctrlds.inp
	echo "dws.txt" >> invctrlds.inp
	echo "50" >> invctrlds.inp
    echo "0.105 1" >> invctrlds.inp
    echo "-0.06 0.06" >> invctrlds.inp

    ./invctrlds



# Relocation

for((k=1;k<=1;k++))
do
# Forward
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

# Prepare input for inverse
	echo "fd.dat" > fd2csm_row.inp
	echo "metafd_row.dat" >> fd2csm_row.inp
	echo "fdcsm_row.dat" >> fd2csm_row.inp
	echo "$ndata" >> fd2csm_row.inp
	./fd2csm_row

	awk '{if(NF==9)print $1,$9}' syndataloc.txt > temporigt
	paste t.txt temporigt | awk '{print ($1+$6)-$5}' > dt2.txt

	echo " Prepare input for inverse finished!"


# Run invloc_svd_wei
    echo "metafd_row.dat" > invloc_svd_wei.inp
    echo "fdcsm_row.dat" >> invloc_svd_wei.inp
    echo "dt2.txt" >> invloc_svd_wei.inp
    echo "out_loc.txt" >> invloc_svd_wei.inp
    echo "$ndata $nevn $nsl" >> invloc_svd_wei.inp
    echo "100000" >> invloc_svd_wei.inp

    ./invloc_svd_wei
    echo "Inversion part finished!!"

	cp syndataloc.txt syndataloctmp_$k.txt

# Convert results of invsvddt to real locations
	echo "syndataloctmp_$k.txt" > dloc2loc_multi.inp
	echo "out_loc.txt" >> dloc2loc_multi.inp
	echo "syndataloc.txt" >> dloc2loc_multi.inp
	echo "$sta_num $nevn" >> dloc2loc_multi.inp
	echo "0.25 0.08" >> dloc2loc_multi.inp
	echo "$minlon $maxlon" >> dloc2loc_multi.inp
	echo "$minlan $maxlan" >> dloc2loc_multi.inp
	echo "$minz $maxz" >> dloc2loc_multi.inp
	echo "$topo_file" >> dloc2loc_multi.inp
	echo "$topo_minlon $topo_maxlon" >> dloc2loc_multi.inp
	echo "$topo_minlan $topo_maxlan" >> dloc2loc_multi.inp
	echo "$topo_xnum $topo_ynum" >> dloc2loc_multi.inp

	./dloc2loc_multi
done
    cp syndataloc.txt syndataloc_$i.txt
	cp out_loc.txt out_loc$i.txt

done

rm paths.dat errt.txt

cp vel.txt ../outf/vel_final.txt
cp syndataloc.txt ../outf/dataloc_out.txt

# exit work dir
cd ../../
pwd
