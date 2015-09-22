#!/bin/bash

reff=initvel_3km1d.txt
inpf1=truevel.txt
inpf2=vel_tp.txt
inpf3=vel_notp.txt
minlon=204.50
maxlon=205.20
minlan=19.41
maxlan=19.41
mindep=0.00
maxdep=25.00
thrshd=0.4


# Files for plotting
num=$(head -1 $reff | awk '{print $1}')
head -2 $inpf2 > truevelper.txt
paste $inpf1 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> truevelper.txt

head -2 $inpf2 > velper_tp.txt
paste $inpf2 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> velper_tp.txt

head -2 $inpf3 > velper_notp.txt
paste $inpf3 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> velper_notp.txt

head -2 $inpf2 > velperdiff_tp.txt
paste $inpf2 $inpf1 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$12}' >> velperdiff_tp.txt

head -2 $inpf3 > velperdiff_notp.txt
paste $inpf3 $inpf1 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$12}' >> velperdiff_notp.txt


# cut slice of truevel
echo "slice_truevel.txt" > cutvdiff_reg1_xyz.inp
echo "truevelper.txt" >> cutvdiff_reg1_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg1_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg1_xyz.inp
echo "$mindep $maxdep" >> cutvdiff_reg1_xyz.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg1_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg1_xyz.inp
echo "203.80 205.383" >> cutvdiff_reg1_xyz.inp
echo "18.80 20.383" >> cutvdiff_reg1_xyz.inp
echo "96 96" >> cutvdiff_reg1_xyz.inp
echo "9.673" >> cutvdiff_reg1_xyz.inp
echo "2.0" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of vel_tp
echo "slice_veltp.txt" > cutvdiff_reg1_xyz.inp
echo "velper_tp.txt" >> cutvdiff_reg1_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg1_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg1_xyz.inp
echo "$mindep $maxdep" >> cutvdiff_reg1_xyz.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg1_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg1_xyz.inp
echo "203.80 205.383" >> cutvdiff_reg1_xyz.inp
echo "18.80 20.383" >> cutvdiff_reg1_xyz.inp
echo "96 96" >> cutvdiff_reg1_xyz.inp
echo "9.673" >> cutvdiff_reg1_xyz.inp
echo "2.0" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of vel_notp
echo "slice_velnotp.txt" > cutvdiff_reg1_xyz.inp
echo "velper_notp.txt" >> cutvdiff_reg1_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg1_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg1_xyz.inp
echo "$mindep $maxdep" >> cutvdiff_reg1_xyz.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg1_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg1_xyz.inp
echo "203.80 205.383" >> cutvdiff_reg1_xyz.inp
echo "18.80 20.383" >> cutvdiff_reg1_xyz.inp
echo "96 96" >> cutvdiff_reg1_xyz.inp
echo "9.673" >> cutvdiff_reg1_xyz.inp
echo "2.0" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of velperdiff_tp
echo "slice_veldifftp.txt" > cutvdiff_reg2_xyz.inp
echo "velperdiff_tp.txt" >> cutvdiff_reg2_xyz.inp
echo "mdres.txt $thrshd" >> cutvdiff_reg2_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg2_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg2_xyz.inp
echo "$mindep $maxdep" >> cutvdiff_reg2_xyz.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg2_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg2_xyz.inp
echo "203.80 205.383" >> cutvdiff_reg2_xyz.inp
echo "18.80 20.383" >> cutvdiff_reg2_xyz.inp
echo "96 96" >> cutvdiff_reg2_xyz.inp
echo "9.673" >> cutvdiff_reg2_xyz.inp
echo "2.0" >> cutvdiff_reg2_xyz.inp

./cutvdiff_reg2_xyz

# cut slice of velperdiff_notp
echo "slice_veldiffnotp.txt" > cutvdiff_reg2_xyz.inp
echo "velperdiff_notp.txt" >> cutvdiff_reg2_xyz.inp
echo "mdres.txt $thrshd" >> cutvdiff_reg2_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg2_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg2_xyz.inp
echo "$mindep $maxdep" >> cutvdiff_reg2_xyz.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg2_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg2_xyz.inp
echo "203.80 205.383" >> cutvdiff_reg2_xyz.inp
echo "18.80 20.383" >> cutvdiff_reg2_xyz.inp
echo "96 96" >> cutvdiff_reg2_xyz.inp
echo "9.673" >> cutvdiff_reg2_xyz.inp
echo "2.0" >> cutvdiff_reg2_xyz.inp

./cutvdiff_reg2_xyz


# cut model resolution slice
echo "slice_mdres.txt" > cut_ressl_reg2_xyz.inp
echo "mdres.txt" >> cut_ressl_reg2_xyz.inp
echo "$reff" >> cut_ressl_reg2_xyz.inp
echo "$minlon $maxlon" >> cut_ressl_reg2_xyz.inp
echo "$minlan $maxlan" >> cut_ressl_reg2_xyz.inp
echo "$mindep $maxdep" >> cut_ressl_reg2_xyz.inp
echo "0.005 0.005 0.5" >> cut_ressl_reg2_xyz.inp
echo "hawaii_tp4tomo2.xyz" >> cut_ressl_reg2_xyz.inp
echo "203.80 205.383" >> cut_ressl_reg2_xyz.inp
echo "18.80 20.383" >> cut_ressl_reg2_xyz.inp
echo "96 96" >> cut_ressl_reg2_xyz.inp
echo "9.673" >> cut_ressl_reg2_xyz.inp
echo "2.0" >> cut_ressl_reg2_xyz.inp

./cut_ressl_reg2_xyz

echo "0.00000 0.00000 0.000000" > elev_poly.txt
cat elev_line.txt >> elev_poly.txt
tail -1 elev_line.txt | awk '{print $1,0.000,0.000}' >> elev_poly.txt
echo "0.00000 0.00000 0.000000" >> elev_poly.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz.bash $thrshd
