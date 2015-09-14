#!/bin/bash


inpf=$1

nmodel=$(head -1 truevel.txt | awk '{print $1}')

paste truevel.txt $inpf | tail -"$nmodel" > tempdiff

# Residual for total nodes
head -2 $inpf > mddiff.txt
awk '{print $1,$2,$3,$8-$4}' tempdiff >>mddiff.txt
awk '{print $8-$4,($8-$4)/$4}' tempdiff > modelresi.txt

rms=$(awk '{a+=($1*$1)}END{print sqrt(a/"'$nmodel'")}' modelresi.txt)
rmsper=$(awk '{a+=($2*$2)}END{print sqrt(a/"'$nmodel'")}' modelresi.txt)

echo $rms $rmsper $inpf "rms and rmsper of diff"
echo $rms $rmsper $inpf "rms and rmsper of diff" >> mdresi.txt

# Residual for nodes with Large DWS value

paste tempdiff  dws.txt | awk '{if($9>100) print $8-$4,($8-$4)/$4}' > modelresi2.txt

nmodel2=$(wc -l modelresi2.txt | awk '{print $1}')
rmsgd=$(awk '{a+=($1*$1)}END{print sqrt(a/"'$nmodel2'")}' modelresi2.txt)
rmsgdper=$(awk '{a+=($2*$2)}END{print sqrt(a/"'$nmodel2'")}' modelresi2.txt)

echo $rmsgd $rmsgdper $inpf "rms and rmsper of diff"
echo $rmsgd $rmsgdper $inpf "rms and rmsper of diff" >> mdresi.txt
