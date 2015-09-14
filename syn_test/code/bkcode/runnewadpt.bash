#!/bin/bash

srcs=35
dlon=0.005
dlan=0.005
dz=0.5
vconst=5.0
vair=1.5
topobtm=9.673
thr0=0.25
thr1=4
thr2=10
dthrshd=0.02
velmd=vz.1D_inverted

# Generate rectangle nodes and tetrahedral grid by tetgen.
velnum=$(head -1 initvel.txt | awk '{print $1}')
echo "$velnum 3 0 0" > step0.node
tail -"$velnum" initvel.txt | awk '{print $1,$2,$3/100}' | cat -n >> step0.node
tetgen -n step0.node
tetgen -rNEFV step0.1.ele > tetstat_step0.txt

echo "3" > node.txt
echo "$velnum" >> node.txt
tail -"$velnum" initvel.txt | awk '{print $1,$2,$3/100}' >> node.txt

tetnum=$(head -1 step0.1.ele | awk '{print $1}')
tetnum2=$( echo " $tetnum + 1 " | bc)
echo "$tetnum" > tethtg.txt
head -$tetnum2 step0.1.ele | tail -$tetnum | awk '{print $2-1,$3-1,$4-1,$5-1}' >>tethtg.txt

echo "$tetnum" > tethnbr.txt
head -$tetnum2 step0.1.neigh | tail -$tetnum | awk '{print $2,$3,$4,$5}' >>tethnbr.txt

# Delete 0 volumn tetrahedrons from tetgen result.

echo "node.txt" > del0vol.inp
echo "tethtg.txt" >> del0vol.inp
echo "teth.txt" >> del0vol.inp
echo "tethvol.txt" >> del0vol.inp
./del0vol


# Run fmm_adpt3d
echo "syndata.txt" > fmm_adpt3d.inp
echo "initvel.txt" >> fmm_adpt3d.inp
echo "$velmd" >> fmm_adpt3d.inp
echo "teth.txt" >> fmm_adpt3d.inp
echo "tethnbr.txt" >> fmm_adpt3d.inp
echo "$srcs" >> fmm_adpt3d.inp
echo "$dlon $dlan $dz" >> fmm_adpt3d.inp
echo "5.0 5.0 5.0" >> fmm_adpt3d.inp
echo "0.045 0.045 5" >> fmm_adpt3d.inp
echo "204.50 205.20" >> fmm_adpt3d.inp
echo "19.20 19.70" >> fmm_adpt3d.inp
echo "0.00 55.00" >> fmm_adpt3d.inp
echo "$thr0 $thr1 $thr2" >> fmm_adpt3d.inp
echo "$dthrshd" >> fmm_adpt3d.inp
echo "2" >> fmm_adpt3d.inp
echo "2" >> fmm_adpt3d.inp
echo "2" >> fmm_adpt3d.inp
echo "$vconst" >> fmm_adpt3d.inp
echo "hawaii_tp4tomo2.xyz" >> fmm_adpt3d.inp
echo "203.80 205.383" >> fmm_adpt3d.inp
echo "18.80 20.383" >> fmm_adpt3d.inp
echo "96 96" >> fmm_adpt3d.inp
echo "$topobtm" >> fmm_adpt3d.inp
echo "$vair" >> fmm_adpt3d.inp

./fmm_adpt3d

cat edge_delgrd.txt newnode.txt > tempvel
wc -l tempvel | awk '{print $1}' > vel.txt
wc -l edge_delgrd.txt | awk '{print $1}' >> vel.txt
cat tempvel >> vel.txt

cp vel.txt vel_0.txt
cp G_1st_norm.txt G_1stnorm_0.txt
cp edge_delgrd.txt edge_delgrd_0.txt
cp newnode.txt newnode_0.txt
cp delgrid.txt delgrid_0.txt

j=0
for((it=1;it<5;it++))
do
	j=$(( $j + 1 ))	
	# Generate rectangle nodes and tetrahedral grid by tetgen.
	velnum=$(head -1 vel.txt | awk '{print $1}')
	echo "$velnum 3 0 0" > step$it.node
	tail -"$velnum" vel.txt | awk '{print $1,$2,$3/100}' | cat -n >> step$it.node
	tetgen -n step$it.node
	tetgen -rNEFV step$it.1.ele > tetstat_step$it.txt

	echo "3" > node.txt
	echo "$velnum" >> node.txt
	tail -"$velnum" vel.txt | awk '{print $1,$2,$3/100}' >> node.txt

	tetnum=$(head -1 step$it.1.ele | awk '{print $1}')
	tetnum2=$( echo " $tetnum + 1 " | bc)
	echo "$tetnum" > tethtg.txt
	head -$tetnum2 step$it.1.ele | tail -$tetnum | awk '{print $2-1,$3-1,$4-1,$5-1}' >>tethtg.txt

	echo "$tetnum" > tethnbr.txt
	head -$tetnum2 step$it.1.neigh | tail -$tetnum | awk '{print $2,$3,$4,$5}' >>tethnbr.txt


# Delete 0 volumn tetrahedrons from tetgen result.

	echo "node.txt" > del0vol.inp
	echo "tethtg.txt" >> del0vol.inp
	echo "teth.txt" >> del0vol.inp
	echo "tethvol.txt" >> del0vol.inp
	./del0vol


# Run fmm_adpt3d
	echo "syndata.txt" > fmm_adpt3d.inp
	echo "vel.txt" >> fmm_adpt3d.inp
	echo "$velmd" >> fmm_adpt3d.inp
	echo "teth.txt" >> fmm_adpt3d.inp
	echo "tethnbr.txt" >> fmm_adpt3d.inp
	echo "$srcs" >> fmm_adpt3d.inp
	echo "$dlon $dlan $dz" >> fmm_adpt3d.inp
	echo "5.0 5.0 5.0" >> fmm_adpt3d.inp
	echo "0.045 0.045 5" >> fmm_adpt3d.inp
	echo "204.50 205.20" >> fmm_adpt3d.inp
	echo "19.20 19.70" >> fmm_adpt3d.inp
	echo "0.00 55.00" >> fmm_adpt3d.inp
	echo "$thr0 $thr1 $thr2" >> fmm_adpt3d.inp
	echo "$dthrshd" >> fmm_adpt3d.inp
	echo "2" >> fmm_adpt3d.inp
	echo "2" >> fmm_adpt3d.inp
	echo "2" >> fmm_adpt3d.inp
	echo "$vconst" >> fmm_adpt3d.inp
	echo "hawaii_tp4tomo2.xyz" >> fmm_adpt3d.inp
	echo "203.80 205.383" >> fmm_adpt3d.inp
	echo "18.80 20.383" >> fmm_adpt3d.inp
	echo "96 96" >> fmm_adpt3d.inp
	echo "$topobtm" >> fmm_adpt3d.inp
	echo "$vair" >> fmm_adpt3d.inp

	./fmm_adpt3d

	cat edge_delgrd.txt newnode.txt > tempvel
	wc -l tempvel | awk '{print $1}' > vel.txt
	wc -l edge_delgrd.txt | awk '{print $1}' >> vel.txt
	cat tempvel >> vel.txt

	cp vel.txt vel_$it.txt
	cp G_1st_norm.txt G_1stnorm_$it.txt
	cp edge_delgrd.txt edge_delgrd_$it.txt
	cp newnode.txt newnode_$it.txt
	cp delgrid.txt delgrid_$it.txt

done
