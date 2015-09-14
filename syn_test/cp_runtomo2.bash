#!/bin/bash


for gsize in 5km
do
	dir=$(echo tomo$gsize)
	echo $dir

	mkdir $dir
	
	mkdir $dir/withouttopo


	# Copy files to tomography directory
	cp ./code/fmm_fw_regul2 ./$dir/withouttopo/
	cp ./code/fmm_fw_regul2.inp ./$dir/withouttopo/
	cp ./code/invsvddt ./$dir/withouttopo/
	cp ./code/invsvddt.inp ./$dir/withouttopo/
	cp ./code/fd2csm ./$dir/withouttopo/
	cp ./code/fd2csm.inp ./$dir/withouttopo/
	cp ./code/ss_regul_topo_tomo2.bash ./$dir/withouttopo/

	cp ./inpf/* ./$dir/withouttopo/

	cp ./gen_grid/vel_$gsize.txt ./$dir/withouttopo/vel.txt

	# Run tomography program
	cd ./$dir/withouttopo/
	pwd
	bash ss_regul_topo_tomo2.bash
	cd ../../
	pwd

	# Make plots
	#mkdir plots/$dir
	#cp ./$dir/withouttopo/out_svd_ds.txt ./plots/$dir/
done

#bash runplots.bash
