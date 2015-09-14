#!/bin/bash


for gsize in 5km
do
	dir=$(echo tomo$gsize)
	echo $dir


	for subdir in withtopo withouttopo
	do
		# Copy files to tomography directory
		cp ./code/fmm_res ./$dir/$subdir/
		cp ./code/ss_regul_topo_res2.bash ./$dir/$subdir/



		# Run tomography program
		cd ./$dir/$subdir/
		pwd
		bash ss_regul_topo_res2.bash
		cd ../../
		pwd
	done

done

