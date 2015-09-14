#!/bin/bash

vfile=$1
nlp=$2
srcs=35
dlon=0.005
dlan=0.005
dz=0.5
vconst=5.0
vair=1.5
topobtm=9.673
thr1=30
datafile=syndataloc.txt
dpmax=25.0
damp=100


# Run fmm_fw_regul2tp
echo "$datafile" > fmm_fw_regul2tp.inp
echo "$vfile" >> fmm_fw_regul2tp.inp
echo "$srcs" >> fmm_fw_regul2tp.inp
echo "$dlon $dlan $dz" >> fmm_fw_regul2tp.inp
echo "5.0 5.0 5.0" >> fmm_fw_regul2tp.inp
echo "0.045 0.045 5" >> fmm_fw_regul2tp.inp
echo "204.50 205.20" >> fmm_fw_regul2tp.inp
echo "19.20 19.70" >> fmm_fw_regul2tp.inp
echo "0.00 $dpmax" >> fmm_fw_regul2tp.inp
echo "2" >> fmm_fw_regul2tp.inp
echo "2" >> fmm_fw_regul2tp.inp
echo "hawaii_tp4tomo2.xyz" >> fmm_fw_regul2tp.inp
echo "203.80 205.383" >> fmm_fw_regul2tp.inp
echo "18.80 20.383" >> fmm_fw_regul2tp.inp
echo "96 96" >> fmm_fw_regul2tp.inp
echo "$topobtm" >> fmm_fw_regul2tp.inp
echo "$vair" >> fmm_fw_regul2tp.inp

./fmm_fw_regul2tp

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

# Run invcgls

echo "locfd.dat" > invcgls.inp
echo "datat.txt" >> invcgls.inp
echo "invlog.txt" >> invcgls.inp
echo "out_loc.txt" >> invcgls.inp
echo "$idata $imodel" >> invcgls.inp
echo "$damp" >> invcgls.inp
echo "1 1" >> invcgls.inp

./invcgls
echo "Inversion part finished!!"


# Run result ctrl

#nlp2=$(echo "$nlp + 1" | bc)
cp vel.txt vel_$nlp.txt
echo "$vfile" > invctrl.inp
echo "out_svd_ds.txt" >> invctrl.inp
echo "vel_$nlp.txt" >> invctrl.inp
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

bash modelresi.bash vel_$nlp.txt
