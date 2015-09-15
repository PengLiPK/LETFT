#!/bin/bash

cp syndataloc_2km10_1dinitcheck.txt syndataloc.txt
cp initvel_3km1d.txt vel.txt
vmodel=vel.txt
damp=100
nevn=2493
nsl=3264
ndata=25465
vair=2.0

cp syndataloc.txt syndataloc_0.txt
for((i=1;i<=20;i++))
do

# Forward
	echo "syndataloc.txt" > fmm_fw_regul2tprl.inp
	echo "$vmodel" >> fmm_fw_regul2tprl.inp
	echo "35" >> fmm_fw_regul2tprl.inp
	echo "0.005 0.005 0.5" >> fmm_fw_regul2tprl.inp
	echo "5.0 5.0 5.0" >> fmm_fw_regul2tprl.inp
	echo "0.045 0.045 5" >> fmm_fw_regul2tprl.inp
	echo "204.50 205.20" >> fmm_fw_regul2tprl.inp
	echo "19.20 19.70" >> fmm_fw_regul2tprl.inp
	echo "0.00 25.00" >> fmm_fw_regul2tprl.inp
	echo "2" >> fmm_fw_regul2tprl.inp
	echo "2" >> fmm_fw_regul2tprl.inp
	echo "hawaii_tp4tomo2.xyz" >> fmm_fw_regul2tprl.inp
	echo "203.80 205.383" >> fmm_fw_regul2tprl.inp
	echo "18.80 20.383" >> fmm_fw_regul2tprl.inp
	echo "96 96" >> fmm_fw_regul2tprl.inp
	echo "9.673" >> fmm_fw_regul2tprl.inp
	echo "$vair" >> fmm_fw_regul2tprl.inp

	./fmm_fw_regul2tprl
	echo "Forward finished!"

# Separate slowness derivative matrix form fd

	awk '{if(NF==9)print $1,$9}' syndataloc.txt > tempt
	paste t.txt tempt | awk '{print ($1+$6)-($5+0.1)}' > dt.txt
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
	echo "$vmodel" >> fmm_fw_regul2tprl.inp
	echo "35" >> fmm_fw_regul2tprl.inp
	echo "0.005 0.005 0.5" >> fmm_fw_regul2tprl.inp
	echo "5.0 5.0 5.0" >> fmm_fw_regul2tprl.inp
	echo "0.045 0.045 5" >> fmm_fw_regul2tprl.inp
	echo "204.50 205.20" >> fmm_fw_regul2tprl.inp
	echo "19.20 19.70" >> fmm_fw_regul2tprl.inp
	echo "0.00 25.00" >> fmm_fw_regul2tprl.inp
	echo "2" >> fmm_fw_regul2tprl.inp
	echo "2" >> fmm_fw_regul2tprl.inp
	echo "hawaii_tp4tomo2.xyz" >> fmm_fw_regul2tprl.inp
	echo "203.80 205.383" >> fmm_fw_regul2tprl.inp
	echo "18.80 20.383" >> fmm_fw_regul2tprl.inp
	echo "96 96" >> fmm_fw_regul2tprl.inp
	echo "9.673" >> fmm_fw_regul2tprl.inp
	echo "$vair" >> fmm_fw_regul2tprl.inp

	./fmm_fw_regul2tprl
	echo "Forward finished!"


# Prepare input for inverse

	echo "fd.dat" > fd2csm_row.inp
	echo "metafd_row.dat" >> fd2csm_row.inp
	echo "fdcsm_row.dat" >> fd2csm_row.inp
	echo "$ndata" >> fd2csm_row.inp
	./fd2csm_row

	awk '{if(NF==9)print $1,$9}' syndataloc.txt > temporigt
	paste t.txt temporigt | awk '{print ($1+$6)-($5+0.1)}' > dt2.txt

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
	echo "35 $nevn" >> dloc2loc_multi.inp
	echo "0.25 0.08" >> dloc2loc_multi.inp
	echo "204.50 205.20" >> dloc2loc_multi.inp
	echo "19.20 19.70" >> dloc2loc_multi.inp
	echo "0.00 25.00" >> dloc2loc_multi.inp
	echo "hawaii_tp4tomo2.xyz" >> dloc2loc_multi.inp
	echo "203.80 205.383" >> dloc2loc_multi.inp
	echo "18.80 20.383" >> dloc2loc_multi.inp
	echo "96 96" >> dloc2loc_multi.inp

	./dloc2loc_multi
done
    cp syndataloc.txt syndataloc_$i.txt
	cp out_loc.txt out_loc$i.txt

done

# Model resolution matrix
echo "dtsep.txt" > mdresolmat.inp
echo "metafd.dat" >> mdresolmat.inp
echo "fdcsm.dat" >> mdresolmat.inp
echo "mdres.txt" >> mdresolmat.inp
echo "$nsldata $nsl" >> mdresolmat.inp
echo "$damp" >> mdresolmat.inp

./mdresolmat
echo "Model resolution maxtrix calculation finished!!"

# Residual
echo "syndataloc.txt" > rmm_fw_regul2tprl.inp
echo "$vmodel" >> rmm_fw_regul2tprl.inp
echo "35" >> rmm_fw_regul2tprl.inp
echo "0.005 0.005 0.5" >> rmm_fw_regul2tprl.inp
echo "5.0 5.0 5.0" >> rmm_fw_regul2tprl.inp
echo "0.045 0.045 5" >> rmm_fw_regul2tprl.inp
echo "204.50 205.20" >> rmm_fw_regul2tprl.inp
echo "19.20 19.70" >> rmm_fw_regul2tprl.inp
echo "0.00 25.00" >> rmm_fw_regul2tprl.inp
echo "2" >> rmm_fw_regul2tprl.inp
echo "2" >> rmm_fw_regul2tprl.inp
echo "hawaii_tp4tomo2.xyz" >> rmm_fw_regul2tprl.inp
echo "203.80 205.383" >> rmm_fw_regul2tprl.inp
echo "18.80 20.383" >> rmm_fw_regul2tprl.inp
echo "96 96" >> rmm_fw_regul2tprl.inp
echo "9.673" >> rmm_fw_regul2tprl.inp
echo "$vair" >> rmm_fw_regul2tprl.inp

./rmm_fw_regul2tprl
echo "Forward finished!"


j=$(( $i - 1 ))
awk '{if(NF==9)print $1,$9}' syndataloc.txt > tempt
paste t.txt tempt | awk '{print ($1+$6)-$5}' > res_$j.txt


rm paths.dat errt.txt
