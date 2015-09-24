#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=syn_work_ps
inpfdir=syn_outf
inpfdir_1dv=input_file
plotdir=syn_plot
#######################################################

# Parameters of input files
#######################################################
reff_p=synvelp_init.txt
inpf1_p=synvelp.txt
inpf2_p=vel_p_final.txt
modresf_p=mdres_p.txt
onedvfile_p=vz.1D_inverted
reff_s=synvels_init.txt
inpf1_s=synvels.txt
inpf2_s=vel_s_final.txt
modresf_s=mdres_s.txt
onedvfile_s=vz.1D_inverted_s
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
topo_vair_p=2.0
topo_vair_s=1.15
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
cp ./$wkdir/$inpfdir/$reff_p ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf1_p ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf2_p ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$modresf_p ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$onedvfile_p ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$reff_s ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf1_s ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$inpf2_s ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$modresf_s ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$onedvfile_s ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$topo_file ./$wkdir/$plotdir/

cp ./src_plot/cutvdiff_reg1_xyz ./$wkdir/$plotdir/
cp ./src_plot/cutvdiff_reg2_xyz ./$wkdir/$plotdir/
cp ./src_plot/cut_ressl_reg2_xyz ./$wkdir/$plotdir/
cp ./src_plot/valmask ./$wkdir/$plotdir/
cp ./src_plot/plot_v_xyz.bash ./$wkdir/$plotdir/


# Enter working dir
cd ./$wkdir/$plotdir/
pwd

# Plotting P velocity model.
# Files for plotting
num=$(head -1 $reff_p | awk '{print $1}')
head -2 $inpf2_p > truevelper.txt
paste $inpf1_p $reff_p | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> truevelper.txt

head -2 $inpf2_p > velper_tp.txt
paste $inpf2_p $reff_p | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> velper_tp.txt

head -2 $inpf2_p > velperdiff_tp.txt
paste $inpf2_p $inpf1_p $reff_p | tail -$num | awk '{print $1,$2,$3,($4-$8)/$12}' >> velperdiff_tp.txt


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
echo "$topo_vair_p" >> cutvdiff_reg1_xyz.inp

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
echo "$topo_vair_p" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of velperdiff_tp
echo "slice_veldifftp.txt" > cutvdiff_reg2_xyz.inp
echo "velperdiff_tp.txt" >> cutvdiff_reg2_xyz.inp
echo "$modresf_p $thrshd" >> cutvdiff_reg2_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg2_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg2_xyz.inp
echo "$minz $maxz" >> cutvdiff_reg2_xyz.inp
echo "$dx $dy $dz" >> cutvdiff_reg2_xyz.inp
echo "$topo_file" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_reg2_xyz.inp
echo "$topo_depth" >> cutvdiff_reg2_xyz.inp
echo "$topo_vair_p" >> cutvdiff_reg2_xyz.inp

./cutvdiff_reg2_xyz

# cut model resolution slice
echo "slice_mdres.txt" > cut_ressl_reg2_xyz.inp
echo "$modresf_p" >> cut_ressl_reg2_xyz.inp
echo "$reff_p" >> cut_ressl_reg2_xyz.inp
echo "$minlon $maxlon" >> cut_ressl_reg2_xyz.inp
echo "$minlan $maxlan" >> cut_ressl_reg2_xyz.inp
echo "$minz $maxz" >> cut_ressl_reg2_xyz.inp
echo "$dx $dy $dz" >> cut_ressl_reg2_xyz.inp
echo "$topo_file" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cut_ressl_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cut_ressl_reg2_xyz.inp
echo "$topo_depth" >> cut_ressl_reg2_xyz.inp
echo "$topo_vair_p" >> cut_ressl_reg2_xyz.inp

./cut_ressl_reg2_xyz

echo "0.00000 0.00000 0.000000" > elev_poly.txt
cat elev_line.txt >> elev_poly.txt
tail -1 elev_line.txt | awk '{print $1,0.000,0.000}' >> elev_poly.txt
echo "0.00000 0.00000 0.000000" >> elev_poly.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz.bash $thrshd $onedvfile_p
mv syn_comp_xyz.ps syn_comp_xyz_p.ps
mv syn_comp_xyz.eps syn_comp_xyz_p.eps


# Plotting S velocity model.
# Files for plotting
num=$(head -1 $reff_s | awk '{print $1}')
head -2 $inpf2_s > truevelper.txt
paste $inpf1_s $reff_s | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> truevelper.txt

head -2 $inpf2_s > velper_tp.txt
paste $inpf2_s $reff_s | tail -$num | awk '{print $1,$2,$3,($4-$8)/$8}' >> velper_tp.txt

head -2 $inpf2_s > velperdiff_tp.txt
paste $inpf2_s $inpf1_s $reff_s | tail -$num | awk '{print $1,$2,$3,($4-$8)/$12}' >> velperdiff_tp.txt


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
echo "$topo_vair_s" >> cutvdiff_reg1_xyz.inp

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
echo "$topo_vair_s" >> cutvdiff_reg1_xyz.inp

./cutvdiff_reg1_xyz


# cut slice of velperdiff_tp
echo "slice_veldifftp.txt" > cutvdiff_reg2_xyz.inp
echo "velperdiff_tp.txt" >> cutvdiff_reg2_xyz.inp
echo "$modresf_s $thrshd" >> cutvdiff_reg2_xyz.inp
echo "$minlon $maxlon" >> cutvdiff_reg2_xyz.inp
echo "$minlan $maxlan" >> cutvdiff_reg2_xyz.inp
echo "$minz $maxz" >> cutvdiff_reg2_xyz.inp
echo "$dx $dy $dz" >> cutvdiff_reg2_xyz.inp
echo "$topo_file" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_reg2_xyz.inp
echo "$topo_depth" >> cutvdiff_reg2_xyz.inp
echo "$topo_vair_s" >> cutvdiff_reg2_xyz.inp

./cutvdiff_reg2_xyz

# cut model resolution slice
echo "slice_mdres.txt" > cut_ressl_reg2_xyz.inp
echo "$modresf_s" >> cut_ressl_reg2_xyz.inp
echo "$reff_s" >> cut_ressl_reg2_xyz.inp
echo "$minlon $maxlon" >> cut_ressl_reg2_xyz.inp
echo "$minlan $maxlan" >> cut_ressl_reg2_xyz.inp
echo "$minz $maxz" >> cut_ressl_reg2_xyz.inp
echo "$dx $dy $dz" >> cut_ressl_reg2_xyz.inp
echo "$topo_file" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cut_ressl_reg2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cut_ressl_reg2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cut_ressl_reg2_xyz.inp
echo "$topo_depth" >> cut_ressl_reg2_xyz.inp
echo "$topo_vair_s" >> cut_ressl_reg2_xyz.inp

./cut_ressl_reg2_xyz

echo "0.00000 0.00000 0.000000" > elev_poly.txt
cat elev_line.txt >> elev_poly.txt
tail -1 elev_line.txt | awk '{print $1,0.000,0.000}' >> elev_poly.txt
echo "0.00000 0.00000 0.000000" >> elev_poly.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz.bash $thrshd $onedvfile_p
mv syn_comp_xyz.ps syn_comp_xyz_s.ps
mv syn_comp_xyz.eps syn_comp_xyz_s.eps

cp *.ps *.eps ../$inpfdir

# Exit dir
cd ../../
