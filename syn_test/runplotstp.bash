#!/bin/bash


for dir in tomo5km
do
	echo $dir
	
	mkdir plots/$dir
	mkdir plots/$dir/withtopo


	# Copy plot files
	cp ./plotcode/*.bash ./plots/$dir/withtopo/
	cp ./plotcode/cutslice_reg2 ./plots/$dir/withtopo/

	cp ./inpf/* ./plots/$dir/withtopo/

	cp ./$dir/withtopo/out_svd_ds.txt ./plots/$dir/withtopo/
	cp ./$dir/withtopo/vel.txt ./plots/$dir/withtopo/

	# Run tomography program
	cd ./plots/$dir/withtopo/
	pwd
	bash run_cut_plot.bash
	cd ../../../
	pwd

done
