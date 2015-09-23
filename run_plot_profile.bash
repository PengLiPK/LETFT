#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=syn_work
inpfdir=syn_outf
inpfdir_1dv=input_file
plotdir=syn_plot
#######################################################

# Parameters of input files
#######################################################
reff=synvel_init.txt
inpf1=synvel.txt
inpf2=vel_final.txt
modresf=mdres.txt
onedvfile=vz.1D_inverted
#######################################################

# Parameters of topography file
#######################################################
topo_file=hawaii_tp4tomo2.xyz
topo_minlon=203.80
topo_maxlon=205.383
topo_minlan=18.80
topo_maxlan=20.383
topo_xnum=96
topo_ynum=96
topo_depth=9.673
topo_vair=2.0
#######################################################

# Parameters of plotting
#######################################################
thrshd=0.4 # Model resolution > 0.4 will be plotted.
minlon=204.50
maxlon=205.20
minlan=19.38
maxlan=19.38
minz=0
maxz=25.0
dx=0.005
dy=0.005
dz=0.5
#######################################################

# Check working dir is exist or not.
if [ ! -d $wkdir ]
then
	echo "There is no dir called $wkdir!"
	echo "You need to run cp_runsyn_vel_t.bash first!"
	exit
fi
if [ ! -d $wkdir/$inpfdir/ ]
then
	echo "There is no dir called $wkdir/$inpfdir/!"
	echo "You need to run cp_runsyn_tomo.bash first!"
	exit
fi
if [ ! -d $wkdir/$plotdir/ ]
then
	mkdir ./$wkdir/$plotdir
fi

# Copy files
cp ./$wkdir/$inpfdir/$reff ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf1 ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf2 ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$modresf ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$onedvfile ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$topo_file ./$wkdir/$plotdir/

cp ./src_plot/cutvdiff_reg1_xyz ./$wkdir/$plotdir/
cp ./src_plot/cutvdiff_reg2_xyz ./$wkdir/$plotdir/
cp ./src_plot/cut_ressl_reg2_xyz ./$wkdir/$plotdir/
cp ./src_plot/valmask ./$wkdir/$plotdir/
cp ./src_plot/plot_v_xyz.bash ./$wkdir/$plotdir/


# Enter working dir
cd ./$wkdir/$plotdir/
pwd

# Files for plotting
num=$(head -1 $reff | awk '{print $1}')
head -2 $inpf2 > truevelper.txt
paste $inpf1 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> truevelper.txt

head -2 $inpf2 > velper_tp.txt
paste $inpf2 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> velper_tp.txt

head -2 $inpf2 > velperdiff_tp.txt
paste $inpf2 $inpf1 $reff | tail -$num | awk '{print $1,$2,$3,($4-$8)/$12}' >> velperdiff_tp.txt


# cut slice of truevel
echo "slice_truevel.txt" > cutvdiff_reg1_xyz.inp
echo "truevelper.txt" >> cutvdiff_reg1_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg1_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg1_xyz.inp
echo "$minz $maxz" >> cutvdiff_reg1_xyz.inp
echo "$dx $dy $dz" >> cutvdiff_reg1_xyz.inp
echo "$topo_file" >> cutvdiff_reg1_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_reg1_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_reg1_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_reg1_xyz.inp
echo "$topo_depth" >> cutvdiff_reg1_xyz.inp
echo "$topo_vair" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of vel_tp
echo "slice_veltp.txt" > cutvdiff_reg1_xyz.inp
echo "velper_tp.txt" >> cutvdiff_reg1_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg1_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg1_xyz.inp
echo "$minz $maxz" >> cutvdiff_reg1_xyz.inp
echo "$dx $dy $dz" >> cutvdiff_reg1_xyz.inp
echo "$topo_file" >> cutvdiff_reg1_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_reg1_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_reg1_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_reg1_xyz.inp
echo "$topo_depth" >> cutvdiff_reg1_xyz.inp
echo "$topo_vair" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of velperdiff_tp
echo "slice_veldifftp.txt" > cutvdiff_reg2_xyz.inp
echo "velperdiff_tp.txt" >> cutvdiff_reg2_xyz.inp
echo "$modresf $thrshd" >> cutvdiff_reg2_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg2_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg2_xyz.inp
echo "$minz $maxz" >> cutvdiff_reg2_xyz.inp
echo "$dx $dy $dz" >> cutvdiff_reg2_xyz.inp
echo "$topo_file" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_reg2_xyz.inp
echo "$topo_depth" >> cutvdiff_reg2_xyz.inp
echo "$topo_vair" >> cutvdiff_reg2_xyz.inp

./cutvdiff_reg2_xyz

# cut model resolution slice
echo "slice_mdres.txt" > cut_ressl_reg2_xyz.inp
echo "$modresf" >> cut_ressl_reg2_xyz.inp
echo "$reff" >> cut_ressl_reg2_xyz.inp
echo "$minlon $maxlon" >> cut_ressl_reg2_xyz.inp
echo "$minlan $maxlan" >> cut_ressl_reg2_xyz.inp
echo "$minz $maxz" >> cut_ressl_reg2_xyz.inp
echo "$dx $dy $dz" >> cut_ressl_reg2_xyz.inp
echo "$topo_file" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cut_ressl_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cut_ressl_reg2_xyz.inp
echo "$topo_depth" >> cut_ressl_reg2_xyz.inp
echo "$topo_vair" >> cut_ressl_reg2_xyz.inp

./cut_ressl_reg2_xyz

echo "0.00000 0.00000 0.000000" > elev_poly.txt
cat elev_line.txt >> elev_poly.txt
tail -1 elev_line.txt | awk '{print $1,0.000,0.000}' >> elev_poly.txt
echo "0.00000 0.00000 0.000000" >> elev_poly.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz.bash $thrshd
cp *.ps *.eps ../$inpfdir

# Exit dir
cd ../../
