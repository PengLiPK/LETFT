#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work
inpfdir=input_file
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
	mkdir $wkdir/syn_tomo/
fi
		

# Copy files to tomography directory
cp ./code/fmm_fw_regul2tp ./$wkdir/syn_tomo/
cp ./code/fmm_fw_regul2tp.inp ./$wkdir/syn_tomo/
cp ./code/invsvddt ./$wkdir/syn_tomo/
cp ./code/invsvddt.inp ./$wkdir/syn_tomo/
cp ./code/fd2csm ./$wkdir/syn_tomo/
cp ./code/fd2csm.inp ./$wkdir/syn_tomo/
cp ./code/ss_regul_topo_tomo2tp.bash ./$wkdir/syn_tomo/

cp ./inpf/* ./$wkdir/syn_tomo/

cp ./gen_grid/vel_$gsize.txt ./$wkdir/syn_tomo/vel.txt

# Run tomography program
cd ./$wkdir/syn_tomo/
pwd
bash ss_regul_topo_tomo2tp.bash
cd ../../
pwd


