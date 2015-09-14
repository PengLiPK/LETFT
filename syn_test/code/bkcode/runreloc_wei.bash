#!/bin/bash

vmodel=vel.txt

damp=0.01

for((i=1;i<=4;i++))
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
	echo "1.5" >> fmm_fw_regul2tprl.inp

	./fmm_fw_regul2tprl
	echo "Forward finished!"


# Prepare input for inverse

	echo "fd.dat" > fd2locfd.inp
	echo "locfd.dat" >> fd2locfd.inp 
	echo "2493" >> fd2locfd.inp
	echo "1" >> fd2locfd.inp
	./fd2locfd

	idata=$(wc -l dt2.txt | awk '{print $1}')

	echo "locfd.dat" > addweight.inp
	echo "locfdw.dat" >> addweight.inp
	echo "dt2.txt" >> addweight.inp
	echo "dt2w.txt" >> addweight.inp
	echo "dt2.txt" >> addweight.inp
	echo "$idata" >> addweight.inp
	echo "1" >> addweight.inp
	./addweight


	echo "locfdw.dat" > fd2csm.inp
	echo "metalocfd.dat" >> fd2csm.inp
	echo "locfdcsm.dat" >> fd2csm.inp
	echo "9972" >> fd2csm.inp
	./fd2csm

	awk '{if(NF==9)print $1}' syndata.txt > tempt
	awk '{if(NF==9)print $9}' syndataloc.txt > temporigt
	paste t.txt tempt temporigt | awk '{print ($1+$6)-($5+0.1)}' > dt2.txt

	echo " Prepare input for inverse finished!"


# Run invcgls
	imodel=$(wc -l metalocfd.dat | awk '{print $1}')

	echo "locfd.dat" > invcgls.inp
	echo "dt2w.txt" >> invcgls.inp
	echo "invlog.txt" >> invcgls.inp
	echo "out_loc.txt" >> invcgls.inp
	echo "$idata $imodel" >> invcgls.inp
	echo "$damp" >> invcgls.inp
	echo "1 1" >> invcgls.inp

	./invcgls
	echo "Inversion part finished!!"


	cp syndataloc.txt syndataloc_$i.txt

# Convert results of invcgls to real locations
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
