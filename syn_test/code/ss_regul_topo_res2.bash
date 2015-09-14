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
nodefile=out_svd_ds.txt
dpmax=55.0

mv t.txt tfwd.txt

head -1 vel.txt > tempds1
imodel=$(head -1 vel.txt | awk '{print $1}')
tail -$imodel vel.txt > tempds2
paste tempds2 $nodefile | awk '{print $1,$2,$3,-$5}' > tempds3
cat tempds1 tempds3 > out_dsl.txt
rm tempds?

# Run fmm_res
echo "$datafile" > fmm_res.inp
echo "out_dsl.txt" >> fmm_res.inp
echo "$oneDvf" >> fmm_res.inp
echo "$srcs" >> fmm_res.inp
echo "$dlon $dlan $dz" >> fmm_res.inp
echo "5.0 5.0 5.0" >> fmm_res.inp
echo "0.045 0.045 5" >> fmm_res.inp
echo "204.50 205.20" >> fmm_res.inp
echo "19.20 19.70" >> fmm_res.inp
echo "0.00 $dpmax" >> fmm_res.inp
echo "2" >> fmm_res.inp
echo "2" >> fmm_res.inp
echo "hawaii_tp4tomo2.xyz" >> fmm_res.inp
echo "203.80 205.383" >> fmm_res.inp
echo "18.80 20.383" >> fmm_res.inp
echo "96 96" >> fmm_res.inp
echo "$topobtm" >> fmm_res.inp
echo "$vair" >> fmm_res.inp

./fmm_res

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

#./fd2csm


# Calculate dt
awk '{if(NF>6)print $0}' $datafile > tempt1
awk '{if(NF>6)print $0}' t.txt > tempt2
paste tempt1 tempt2 | awk '{print $9-$1}' > res.txt

resstd=$(awk '{a+=($1*$1)}END{print sqrt(a/"'$idata'")}' res.txt)
cdir=$(pwd)
echo $resstd $cdir >> ../../resstd.txt


