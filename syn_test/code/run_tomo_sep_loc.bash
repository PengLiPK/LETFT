#!/bin/bash

#cp syndataloc_initcheck.txt syndataloc.txt
#cp initvel.txt vel.txt
vmodel=vel.txt
damp=0
nevn=2493
nsl=990
ndata=25465

for((i=3;i<=5;i++))
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
	echo "3.5" >> fmm_fw_regul2tprl.inp

	./fmm_fw_regul2tprl
	echo "Forward finished!"

# Separate slowness derivative matrix form fd

	awk '{if(NF==9)print $1,$9}' syndataloc.txt > tempt
	paste t.txt tempt | awk '{print ($1+$6)-($5+0.1)}' > dt.txt
	#paste t.txt tempt | awk '{print $1-$5}' > dt.txt
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


# Run invsvddt
	echo "dtsep.txt" > invsvddt.inp
	echo "metafd.dat" >> invsvddt.inp
	echo "fdcsm.dat" >> invsvddt.inp
	echo "out_dsl.txt" >> invsvddt.inp
	echo "$nsldata $nsl" >> invsvddt.inp
		
	./invsvddt
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
    echo "-0.11 0.11" >> invctrlds.inp

    ./invctrlds



# Relocation
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
	echo "3.5" >> fmm_fw_regul2tprl.inp

	./fmm_fw_regul2tprl
	echo "Forward finished!"


# Prepare input for inverse

	echo "fd.dat" > fd2locfd.inp
	echo "locfd.dat" >> fd2locfd.inp 
	echo "$nevn" >> fd2locfd.inp
	echo "1" >> fd2locfd.inp
	./fd2locfd

	para1=$(( $nevn * 4 ))
	echo "locfd.dat" > fd2csm.inp
	echo "metalocfd.dat" >> fd2csm.inp
	echo "locfdcsm.dat" >> fd2csm.inp
	echo "$para1" >> fd2csm.inp
	./fd2csm

	awk '{if(NF==9)print $1}' syndata.txt > tempt
	awk '{if(NF==9)print $9}' syndataloc.txt > temporigt
	paste t.txt tempt temporigt | awk '{print ($1+$6)-($5+0.1)}' > dt2.txt

	echo " Prepare input for inverse finished!"


# Run invcgls
	idata=$(wc -l dt2.txt | awk '{print $1}')
	imodel=$(wc -l metalocfd.dat | awk '{print $1}')

	echo "locfd.dat" > invcgls.inp
	echo "dt2.txt" >> invcgls.inp
	echo "invlog.txt" >> invcgls.inp
	echo "out_loc.txt" >> invcgls.inp
	echo "$idata $imodel" >> invcgls.inp
	echo "$damp" >> invcgls.inp
	echo "1 1" >> invcgls.inp

	./invcgls
	echo "Inversion part finished!!"


	cp syndataloc.txt syndataloc_$i.txt
	cp out_loc.txt out_loc$i.txt

# Convert results of invsvddt to real locations
	echo "syndataloc_$i.txt" > dloc2loc_multi.inp
	echo "out_loc.txt" >> dloc2loc_multi.inp
	echo "syndataloc.txt" >> dloc2loc_multi.inp
	echo "35 2493" >> dloc2loc_multi.inp
	echo "1 0.3" >> dloc2loc_multi.inp
	echo "204.50 205.20" >> dloc2loc_multi.inp
	echo "19.20 19.70" >> dloc2loc_multi.inp
	echo "0.00 25.00" >> dloc2loc_multi.inp
	echo "hawaii_tp4tomo2.xyz" >> dloc2loc_multi.inp
	echo "203.80 205.383" >> dloc2loc_multi.inp
	echo "18.80 20.383" >> dloc2loc_multi.inp
	echo "96 96" >> dloc2loc_multi.inp

	./dloc2loc_multi

done
