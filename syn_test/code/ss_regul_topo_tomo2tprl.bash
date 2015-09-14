#!/bin/bash

srcs=35
dlon=0.005
dlan=0.005
dz=0.5
vconst=5.0
vair=1.5
topobtm=9.673
thr1=30
datafile=syndata.txt
oneDvf=vz.1D_inverted
nodefile=vel.txt
dpmax=55.0


# Run fmm_fw_regul2tprl
echo "$datafile" > fmm_fw_regul2tprl.inp
echo "$nodefile" >> fmm_fw_regul2tprl.inp
echo "$oneDvf" >> fmm_fw_regul2tprl.inp
echo "$srcs" >> fmm_fw_regul2tprl.inp
echo "$dlon $dlan $dz" >> fmm_fw_regul2tprl.inp
echo "5.0 5.0 5.0" >> fmm_fw_regul2tprl.inp
echo "0.045 0.045 5" >> fmm_fw_regul2tprl.inp
echo "204.50 205.20" >> fmm_fw_regul2tprl.inp
echo "19.20 19.70" >> fmm_fw_regul2tprl.inp
echo "0.00 $dpmax" >> fmm_fw_regul2tprl.inp
echo "2" >> fmm_fw_regul2tprl.inp
echo "2" >> fmm_fw_regul2tprl.inp
echo "hawaii_tp4tomo2.xyz" >> fmm_fw_regul2tprl.inp
echo "203.80 205.383" >> fmm_fw_regul2tprl.inp
echo "18.80 20.383" >> fmm_fw_regul2tprl.inp
echo "96 96" >> fmm_fw_regul2tprl.inp
echo "$topobtm" >> fmm_fw_regul2tprl.inp
echo "$vair" >> fmm_fw_regul2tprl.inp

./fmm_fw_regul2tprl

rm tsource.dat
rm err00*.txt
rm paths.dat

imodel=$(head -1 vel.txt | awk '{print $1}')
idata=$(wc -l t.txt | awk '{print $1}')

# Run fd2csm

echo "fd.dat" > fd2csm.inp
echo "metafd.dat" >> fd2csm.inp
echo "fdcsm.dat" >> fd2csm.inp
echo "$imodel" >> fd2csm.inp

./fd2csm


# Calculate dt
awk '{if(NF>6)print $0}' $datafile > tempt
paste t.txt tempt | awk '{print $1-$5}' > dt.txt


# Run invsvddt

echo "dt.txt" > invsvddt.inp
echo "metafd.dat" >> invsvddt.inp
echo "fdcsm.dat" >> invsvddt.inp
echo "out_svd_ds.txt" >> invsvddt.inp
echo "$idata $imodel" >> invsvddt.inp

./invsvddt
