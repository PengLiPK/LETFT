#!/bin/bash


nmodel=$(head -1 vel.txt | awk '{print $1}')

tail -"$nmodel" vel.txt > tempv


# Residual for total nodes
paste tempv out_svd_ds_ctrl.txt | awk '{print (1/$5)-$4,((1/$5)-$4)/$4}' > modelresi.txt

awk '{a+=($1*$1)}END{print sqrt(a/"'$nmodel'")}' modelresi.txt
awk '{a+=($2*$2)}END{print sqrt(a/"'$nmodel'")}' modelresi.txt

# Residual for nodes with Large DWS value

paste tempv out_svd_ds_ctrl.txt dwstemp.txt | awk '{if($6>100) print (1/$5)-$4,((1/$5)-$4)/$4}' > modelresi2.txt

nmodel2=$(wc -l modelresi2.txt | awk '{print $1}')
awk '{a+=($1*$1)}END{print sqrt(a/"'$nmodel2'")}' modelresi2.txt
awk '{a+=($2*$2)}END{print sqrt(a/"'$nmodel2'")}' modelresi2.txt
