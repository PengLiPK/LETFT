#!/bin/bash


for dir in tomo5km
do
	echo $dir
	
	mkdir plots/$dir
	mkdir plots/$dir/withtopob


	# Copy plot files
	cp ./plotcode/*.bash ./plots/$dir/withtopob/
	cp ./plotcode/cutslice_reg2b ./plots/$dir/withtopob/

	cp ./inpf/* ./plots/$dir/withtopob/

	cp ./$dir/withtopo/out_svd_ds.txt ./plots/$dir/withtopob/
	cp ./$dir/withtopo/vel.txt ./plots/$dir/withtopob/

	# Run tomography program
	cd ./plots/$dir/withtopob/
	pwd
	bash run_cut_plotb.bash
	cd ../../../
	pwd

done
