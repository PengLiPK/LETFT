#!/bin/bash

#for((num=21;num<=21;num++))
#do
#vfile=vel_$num.txt
vfile=vel.txt

for dir in tomo3km_2k_10p1d
do
	echo $dir
	
	mkdir plots/$dir


	# Copy plot files
	cp ./plotcode/*.bash ./plots/$dir/
	cp ./plotcode/cutv_reg2 ./plots/$dir/
	cp ./plotcode/cutvdiff_reg2 ./plots/$dir/
	cp ./plotcode/cut_ressl_reg2 ./plots/$dir/

	cp ./inpf/* ./plots/$dir/

	cp ./$dir/$vfile ./plots/$dir/
	cp ./$dir/mddiff.txt ./plots/$dir/
	cp ./$dir/mddiffper.txt ./plots/$dir/
	cp ./$dir/mdres.txt ./plots/$dir/
	cp ./$dir/dws.txt ./plots/$dir/
	cp ./$dir/rct.txt ./plots/$dir/

	# Run tomography program
	cd ./plots/$dir/
	pwd
	bash run_cutv_plot.bash $vfile

	cp out_svd_n19.4v.txt.ps model2.$vfile.n19.4.ps
	cd ../../
	pwd

done
#done
