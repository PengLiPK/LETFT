#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work
inpfdir=input_file
#######################################################

# Parameters of topography file
#######################################################
topo_file=hawaii_tp4tomo2.xyz
initvel_file=initvel_3km1d.txt
data_file=syndataloc_2km10_1dinitcheck.txt
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
		
# Copy input files
cp ./$inpfdir/$topo_file ./$wkdir/syn_tomo/
cp ./$inpfdir/$initvel_file ./$wkdir/syn_tomo/
cp ./$inpfdir/$data_file ./$wkdir/syn_tomo/


# Copy files to tomography directory
cp ./src/fmm_fw_regul2tprl ./$wkdir/syn_tomo/
cp ./src/fd2csm_row ./$wkdir/syn_tomo/
cp ./src/fdsep_qr ./$wkdir/syn_tomo/
cp ./src/fd2csm ./$wkdir/syn_tomo/
cp ./src/calmdres2 ./$wkdir/syn_tomo/
cp ./src/invsvddt_damp ./$wkdir/syn_tomo/
cp ./src/invctrlds ./$wkdir/syn_tomo/
cp ./src/invloc_svd_wei ./$wkdir/syn_tomo/
cp ./src/dloc2loc_multi ./$wkdir/syn_tomo/
cp ./src/mdresolmat ./$wkdir/syn_tomo/

cp ./run_scripts/run_tomo_sep_loc.bash ./$wkdir/syn_tomo/


# Run tomography program
cd ./$wkdir/syn_tomo/
pwd
bash run_tomo_sep_loc.bash
cd ../../
pwd


