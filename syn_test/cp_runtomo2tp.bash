#!/bin/bash


for gsize in 5km
do
	dir=$(echo tomo$gsize)
	echo $dir

	mkdir $dir
	
	mkdir $dir/withtopo


	# Copy files to tomography directory
	cp ./code/fmm_fw_regul2tp ./$dir/withtopo/
	cp ./code/fmm_fw_regul2tp.inp ./$dir/withtopo/
	cp ./code/invsvddt ./$dir/withtopo/
	cp ./code/invsvddt.inp ./$dir/withtopo/
	cp ./code/fd2csm ./$dir/withtopo/
	cp ./code/fd2csm.inp ./$dir/withtopo/
	cp ./code/ss_regul_topo_tomo2tp.bash ./$dir/withtopo/

	cp ./inpf/* ./$dir/withtopo/

	cp ./gen_grid/vel_$gsize.txt ./$dir/withtopo/vel.txt

	# Run tomography program
	cd ./$dir/withtopo/
	pwd
	bash ss_regul_topo_tomo2tp.bash
	cd ../../
	pwd

	# Make plots
	#mkdir plots/$dir
	#cp ./$dir/withtopo/out_svd_ds.txt ./plots/$dir/
done

#bash runplots.bash
