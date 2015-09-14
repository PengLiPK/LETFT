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
dpmax=25.0


# Run initinpv
echo "initvel.txt" > initinpv.inp
echo "vel.txt" >> initinpv.inp
echo "$oneDvf" >> initinpv.inp

./initinpv
echo "initinpv finished!!"


# Run fmm_fw_regul2
echo "$datafile" > fmm_fw_regul2.inp
echo "initvel.txt" >> fmm_fw_regul2.inp
echo "$srcs" >> fmm_fw_regul2.inp
echo "$dlon $dlan $dz" >> fmm_fw_regul2.inp
echo "5.0 5.0 5.0" >> fmm_fw_regul2.inp
echo "0.045 0.045 5" >> fmm_fw_regul2.inp
echo "204.50 205.20" >> fmm_fw_regul2.inp
echo "19.20 19.70" >> fmm_fw_regul2.inp
echo "0.00 $dpmax" >> fmm_fw_regul2.inp
echo "2" >> fmm_fw_regul2.inp
echo "2" >> fmm_fw_regul2.inp
echo "hawaii_tp4tomo2.xyz" >> fmm_fw_regul2.inp
echo "203.80 205.383" >> fmm_fw_regul2.inp
echo "18.80 20.383" >> fmm_fw_regul2.inp
echo "96 96" >> fmm_fw_regul2.inp
echo "$topobtm" >> fmm_fw_regul2.inp
echo "$vair" >> fmm_fw_regul2.inp

./fmm_fw_regul2

rm tsource.dat
rm err*.txt
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
cp tempt datat.txt

# Run invsvddt

echo "datat.txt" > invsvddt.inp
echo "metafd.dat" >> invsvddt.inp
echo "fdcsm.dat" >> invsvddt.inp
echo "out_svd_ds.txt" >> invsvddt.inp
echo "$idata $imodel" >> invsvddt.inp

./invsvddt
echo "Inversion part finished!!"

# Run result ctrl
echo "vel_1.txt" > invctrl.inp
echo "out_svd_ds.txt" >> invctrl.inp
echo "initvel.txt" >> invctrl.inp
echo "0.105 1" >> invctrl.inp
echo "-0.1 0.1" >> invctrl.inp

./invctrl
echo "Results ctrl finished!!"


# Run calmdres2, calculating DWS values, ray counts
echo "metafd.dat" > calmdres2.inp
echo "fdcsm.dat" >> calmdres2.inp
echo "rct.txt" >> calmdres2.inp
echo "dws.txt" >> calmdres2.inp
echo "$idata $imodel" >> calmdres2.inp
    
./calmdres2
echo "Run calmdres2 finished!!"

bash modelresi.bash vel_1.txt
