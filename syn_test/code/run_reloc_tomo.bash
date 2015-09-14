#!/bin/bash


cp initvel.txt vel.txt
cp syndataloc_init.txt syndataloc.txt
for((i=1;i<=4;i++))
do
	echo $i start!!
	bash runreloc.bash
	cp syndataloc.txt syndataloc_iter$i.txt
	bash ss_regul_topo_tomo2tploop.bash vel.txt $i
done
