#!/bin/bash


for((i=0;i<=24;i++)) 
do 
	bash tres.bash res_$i.txt; 
done

for((i=0;i<=25;i++)) 
do 
	bash locresi.bash syndataloc_$i.txt; 
done
bash locresi.bash syndataloc.txt; 


for((i=1;i<=25;i++)) 
do 
	bash modelresi.bash vel_$i.txt
done
bash modelresi.bash vel.txt
