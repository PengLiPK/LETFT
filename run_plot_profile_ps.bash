#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work_p
inpfdir=outf
inpfdir_1dv=input_file
plotdir=plot
#######################################################

# Parameters of input files
#######################################################
inpf_p=vel_p_final.txt
modresf_p=mdres_p.txt
inpf_s=vel_s_final.txt
modresf_s=mdres_s.txt
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
minz=0
maxz=25.0
dx=0.5
dy=0.5
lineseg=3
x[1]=204.66; y[1]=19.46
x[2]=204.8; y[2]=19.35
x[3]=205.2; y[3]=19.52
#######################################################

# Check working dir is exist or not.
if [ ! -d $wkdir ]
then
	echo "There is no dir called $wkdir!"
	echo "You need to run cp_run_init.bash first!"
	exit
fi
if [ ! -d $wkdir/$inpfdir/ ]
then
	echo "There is no dir called $wkdir/$inpfdir/!"
	echo "You need to run cp_run_tomo.bash first!"
	exit
fi
if [ ! -d $wkdir/$plotdir/ ]
then
	mkdir ./$wkdir/$plotdir
fi

# Copy files
cp ./$wkdir/$inpfdir/$inpf1 ./$wkdir/$plotdir/
cp ./$wkdir/$inpfdir/$modresf ./$wkdir/$plotdir/
cp ./$inpfdir_1dv/$topo_file ./$wkdir/$plotdir/

cp ./src_plot/cutvdiff_line2_xyz ./$wkdir/$plotdir/
cp ./src_plot/valmask ./$wkdir/$plotdir/
cp ./src_plot/plot_v_xyz_line.bash ./$wkdir/$plotdir/


# Enter working dir
cd ./$wkdir/$plotdir/
pwd

# Files for plotting
echo "$lineseg" > tmpsegs.txt
for((i=1;i<=$lineseg;i++))
do
	echo "${x[$i]} ${y[$i]}" >> tmpsegs.txt
done

# cut slice of truevel
echo "slice_veltp_1.txt" > cutvdiff_line2_xyz.inp
echo "slice_mdres_1.txt" >> cutvdiff_line2_xyz.inp
echo "ndprj_1.txt 1.5" >> cutvdiff_line2_xyz.inp
echo "$inpf_p" >> cutvdiff_line2_xyz.inp
echo "$modresf_p $thrshd" >> cutvdiff_line2_xyz.inp
echo "tmpsegs.txt" >> cutvdiff_line2_xyz.inp
echo "$minz $maxz" >> cutvdiff_line2_xyz.inp
echo "$dx $dy" >> cutvdiff_line2_xyz.inp
echo "$topo_file" >> cutvdiff_line2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_line2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_line2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_line2_xyz.inp
echo "$topo_depth" >> cutvdiff_line2_xyz.inp
echo "$topo_vair_p" >> cutvdiff_line2_xyz.inp
./cutvdiff_line2_xyz
mv elev_line.txt elev_line1.txt
mv elev_poly.txt elev_poly1.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz_line.bash $thrshd
mv veltp_line.ps vel_p_line.ps
mv veltp_line.eps vel_p_line.eps

# cut slice of truevel
echo "slice_veltp_1.txt" > cutvdiff_line2_xyz.inp
echo "slice_mdres_1.txt" >> cutvdiff_line2_xyz.inp
echo "ndprj_1.txt 1.5" >> cutvdiff_line2_xyz.inp
echo "$inpf_s" >> cutvdiff_line2_xyz.inp
echo "$modresf_s $thrshd" >> cutvdiff_line2_xyz.inp
echo "tmpsegs.txt" >> cutvdiff_line2_xyz.inp
echo "$minz $maxz" >> cutvdiff_line2_xyz.inp
echo "$dx $dy" >> cutvdiff_line2_xyz.inp
echo "$topo_file" >> cutvdiff_line2_xyz.inp
echo "$topo_minlon $topo_maxlon" >> cutvdiff_line2_xyz.inp
echo "$topo_minlan $topo_maxlan" >> cutvdiff_line2_xyz.inp
echo "$topo_xnum $topo_ynum" >> cutvdiff_line2_xyz.inp
echo "$topo_depth" >> cutvdiff_line2_xyz.inp
echo "$topo_vair_s" >> cutvdiff_line2_xyz.inp
./cutvdiff_line2_xyz
mv elev_line.txt elev_line1.txt
mv elev_poly.txt elev_poly1.txt

# Plot velocity and pertubation of velocity
bash plot_v_xyz_line.bash $thrshd
mv veltp_line.ps vel_s_line.ps
mv veltp_line.eps vel_s_line.eps

cp *.ps *.eps ../$inpfdir

# Exit dir
cd ../../
