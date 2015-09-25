#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work_ps
inpfdir=input_file
#######################################################

# Parameters of input file
######################################################
initvel_file_p=init_vel_p.txt
data_file_p=p_data.txt
initvel_file_s=init_vel_s.txt
data_file_s=s_data.txt
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
topo_vair_s=1.15
#######################################################

# Parameters of tomography
#######################################################
dampp=200
damps=150
nevn=1800
npmodel=3264
npdata=30102
nsdata=12605
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
	echo "You need to run cp_run_init_ps.bash first!"
	exit
fi
if [ ! -d $wkdir/tomo/ ]
then
	mkdir $wkdir/tomo/
fi
		
# Copy input files
cp ./$inpfdir/$topo_file ./$wkdir/tomo/
cp ./$wkdir/outf/$initvel_file_p ./$wkdir/tomo/
cp ./$wkdir/outf/$data_file_p ./$wkdir/tomo/
cp ./$wkdir/outf/$initvel_file_s ./$wkdir/tomo/
cp ./$wkdir/outf/$data_file_s ./$wkdir/tomo/


# Copy files to tomography directory
cp ./src/fmm_fw_regul2tprl ./$wkdir/tomo/
cp ./src/comb_ps ./$wkdir/tomo/
cp ./src/fd2csm_row ./$wkdir/tomo/
cp ./src/fdsep_qr_wei ./$wkdir/tomo/
cp ./src/fd2csm ./$wkdir/tomo/
cp ./src/calmdres2 ./$wkdir/tomo/
cp ./src/invsvddt_damp_ps ./$wkdir/tomo/
cp ./src/invctrlds ./$wkdir/tomo/
cp ./src/invloc_svd_wei ./$wkdir/tomo/
cp ./src/dloc2loc_multi ./$wkdir/tomo/


# Run tomography program
# Enter work dir
cd ./$wkdir/tomo/
pwd
mv $data_file_p syndataloc_p.txt
mv $initvel_file_p vel_p.txt
mv $data_file_s syndataloc_s.txt
mv $initvel_file_s vel_s.txt

ndata=$(( $npdata + $nsdata ))
nmodel=$(( $npmodel + $npmodel ))

cp syndataloc_p.txt syndataloc_p_0.txt
cp syndataloc_s.txt syndataloc_s_0.txt
for((i=1;i<=2;i++))
do

# Forward
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

# Combile s and p's fd and t into single matrix and vector.
	echo "fd_p.dat t_p.txt" > comb_ps.inp
	echo "fd_s.dat t_s.txt" >> comb_ps.inp
	echo "$npdata $npmodel $nevn" >> comb_ps.inp
	echo "fd.dat t.txt" >> comb_ps.inp
	./comb_ps

# Separate slowness derivative matrix form fd

	cat syndataloc_p.txt syndataloc_s.txt | awk '{if(NF==9)print $1,$9}' > tempt
	paste t.txt tempt | awk '{print ($1+$3)-$2}' > dt.txt
	j=$(( $i - 1 ))
	cp dt.txt res_$j.txt
	resrms2=$(awk '{a+=$1*$1}END{print sqrt(a/"'$ndata'")}' res_$j.txt)
	echo $resrms2 >> statres_rms.txt

	echo "fd.dat" > fd2csm_row.inp
	echo "metafd_row.dat" >> fd2csm_row.inp
	echo "fdcsm_row.dat" >> fd2csm_row.inp
	echo "$ndata" >> fd2csm_row.inp
	./fd2csm_row

	
	echo "metafd_row.dat" > fdsep_qr_wei.inp
	echo "fdcsm_row.dat" >> fdsep_qr_wei.inp
	echo "dt.txt" >> fdsep_qr_wei.inp
	echo "fdloc.dat" >> fdsep_qr_wei.inp
	echo "dtloc.txt" >> fdsep_qr_wei.inp
	echo "rloc.dat" >> fdsep_qr_wei.inp
	echo "fdsep.dat" >> fdsep_qr_wei.inp
	echo "dtsep.txt" >> fdsep_qr_wei.inp
	echo "$ndata $nevn $nmodel" >> fdsep_qr_wei.inp
	./fdsep_qr_wei
   

# Prepare input for inverse

	echo "fdsep.dat" > fd2csm.inp
	echo "metafd.dat" >> fd2csm.inp
	echo "fdcsm.dat" >> fd2csm.inp
	echo "$nmodel" >> fd2csm.inp
	./fd2csm


	nsldata=$(wc -l dtsep.txt | awk '{print $1}')

# Run calmdres2, calculating DWS values, ray counts
    echo "metafd.dat" > calmdres2.inp
    echo "fdcsm.dat" >> calmdres2.inp
    echo "rct.txt" >> calmdres2.inp
    echo "dws.txt" >> calmdres2.inp
    echo "$nsldata $nmodel" >> calmdres2.inp
  
    ./calmdres2
	head -$npmodel dws.txt > dws_p.txt
	tail -$npmodel dws.txt > dws_s.txt
	head -$npmodel rct.txt > rct_p.txt
	tail -$npmodel rct.txt > rct_s.txt
    echo "Run calmdres2 finished!!"

	echo " Prepare input for inverse finished!"


# Run invsvddt_damp_ps
    echo "dtsep.txt" > invsvddt_damp_ps.inp
    echo "metafd.dat" >> invsvddt_damp_ps.inp
	echo "fdcsm.dat" >> invsvddt_damp_ps.inp
    echo "out_dsl.txt" >> invsvddt_damp_ps.inp
    echo "$nsldata $npmodel $npmodel" >> invsvddt_damp_ps.inp
    echo "$dampp $damps" >> invsvddt_damp_ps.inp

    ./invsvddt_damp_ps


	echo "Inversion part finished!!"
	
	cp out_dsl.txt out_dsl$i.txt
	head -$npmodel out_dsl.txt > out_dsl_p.txt
	tail -$npmodel out_dsl.txt > out_dsl_s.txt
	cp vel_p.txt vel_p_$i.txt
	cp vel_s.txt vel_s_$i.txt

# Convert results of invsvddt to real velocity
    echo "vel_p.txt" > invctrlds.inp
	echo "out_dsl_p.txt" >> invctrlds.inp
	echo "vel_p_$i.txt" >> invctrlds.inp
	echo "dws_p.txt" >> invctrlds.inp
	echo "50" >> invctrlds.inp
    echo "0.101 1" >> invctrlds.inp
    echo "-0.06 0.06" >> invctrlds.inp
    ./invctrlds

    echo "vel_s.txt" > invctrlds.inp
	echo "out_dsl_s.txt" >> invctrlds.inp
	echo "vel_s_$i.txt" >> invctrlds.inp
	echo "dws_s.txt" >> invctrlds.inp
	echo "50" >> invctrlds.inp
    echo "0.101 1" >> invctrlds.inp
    echo "-0.06 0.06" >> invctrlds.inp
    ./invctrlds


# Relocation

for((k=1;k<=1;k++))
do
# Forward
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

# Combine s and p fd and t into a single matrix and vector.
	echo "fd_p.dat t_p.txt" > comb_ps.inp
	echo "fd_s.dat t_s.txt" >> comb_ps.inp
	echo "$npdata $npmodel $nevn" >> comb_ps.inp
	echo "fd.dat t.txt" >> comb_ps.inp
	./comb_ps

# Prepare input for inverse
	echo "fd.dat" > fd2csm_row.inp
	echo "metafd_row.dat" >> fd2csm_row.inp
	echo "fdcsm_row.dat" >> fd2csm_row.inp
	echo "$ndata" >> fd2csm_row.inp
	./fd2csm_row

	cat syndataloc_p.txt syndataloc_s.txt |  awk '{if(NF==9)print $1,$9}' > tempt
	paste t.txt tempt | awk '{print ($1+$3)-$2}' > dt2.txt

	echo " Prepare input for inverse finished!"


# Run invloc_svd_wei
    echo "metafd_row.dat" > invloc_svd_wei.inp
    echo "fdcsm_row.dat" >> invloc_svd_wei.inp
    echo "dt2.txt" >> invloc_svd_wei.inp
    echo "out_loc.txt" >> invloc_svd_wei.inp
    echo "$ndata $nevn $nmodel" >> invloc_svd_wei.inp
    echo "100" >> invloc_svd_wei.inp

    ./invloc_svd_wei
    echo "Inversion part finished!!"

	cp syndataloc_p.txt syndataloctmp_p_$k.txt
	cp syndataloc_s.txt syndataloctmp_s_$k.txt

# Convert results of invsvddt to real locations
	echo "syndataloctmp_p_$k.txt" > dloc2loc_multi.inp
	echo "out_loc.txt" >> dloc2loc_multi.inp
	echo "syndataloc_p.txt" >> dloc2loc_multi.inp
	echo "$sta_num_p $nevn" >> dloc2loc_multi.inp
	echo "0.25 0.08" >> dloc2loc_multi.inp
	echo "$minlon $maxlon" >> dloc2loc_multi.inp
	echo "$minlan $maxlan" >> dloc2loc_multi.inp
	echo "$minz $maxz" >> dloc2loc_multi.inp
	echo "$topo_file" >> dloc2loc_multi.inp
	echo "$topo_minlon $topo_maxlon" >> dloc2loc_multi.inp
	echo "$topo_minlan $topo_maxlan" >> dloc2loc_multi.inp
	echo "$topo_xnum $topo_ynum" >> dloc2loc_multi.inp
	./dloc2loc_multi

	echo "syndataloctmp_s_$k.txt" > dloc2loc_multi.inp
	echo "out_loc.txt" >> dloc2loc_multi.inp
	echo "syndataloc_s.txt" >> dloc2loc_multi.inp
	echo "$sta_num_s $nevn" >> dloc2loc_multi.inp
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
    cp syndataloc_p.txt syndataloc_p_$i.txt
    cp syndataloc_s.txt syndataloc_s_$i.txt
	cp out_loc.txt out_loc$i.txt

done

rm paths.dat errt.txt

cp vel_p.txt ../outf/vel_p_final.txt
cp vel_s.txt ../outf/vel_s_final.txt
cp syndataloc_p.txt ../outf/dataloc_p_out.txt
cp syndataloc_s.txt ../outf/dataloc_s_out.txt

# exit work dir
cd ../../
pwd
