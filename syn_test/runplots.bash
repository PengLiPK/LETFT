#!/bin/bash


for dir in tomo5km
do
	echo $dir
	
	mkdir plots/$dir
	mkdir plots/$dir/withouttopo


	# Copy plot files
	cp ./plotcode/*.bash ./plots/$dir/withouttopo/
	cp ./plotcode/cutslice_reg2 ./plots/$dir/withouttopo/

	cp ./inpf/* ./plots/$dir/withouttopo/

	cp ./$dir/withouttopo/out_svd_ds.txt ./plots/$dir/withouttopo/
	cp ./$dir/withouttopo/vel.txt ./plots/$dir/withouttopo/

	# Run tomography program
	cd ./plots/$dir/withouttopo/
	pwd
	bash run_cut_plot.bash
	cd ../../../
	pwd

done
