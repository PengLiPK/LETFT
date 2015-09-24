#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work_p
inpfdir=input_file
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

# Parameters of initial velocity model
#######################################################
node_structure=vstrct_tp3km.txt
oned_vel=vz.1D_inverted
#######################################################

# Parameters of migrate locations of events which are
# above the Earth's surface
#######################################################
data_file=p_data.txt
sta_num=35
#######################################################

# Check working dir is exist or not.
if [ ! -d $wkdir ]
then
	mkdir $wkdir
fi

if [ ! -d $wkdir/init_vel_t/ ]
then
	mkdir $wkdir/init_vel_t/
fi

if [ ! -d $wkdir/outf/ ]
then
	mkdir $wkdir/outf/
fi

# Enter $wkdir/init_vel_t
cd $wkdir/init_vel_t

# Copy input files
cp ../../$inpfdir/$node_structure .
cp ../../$inpfdir/$topo_file .
cp ../../$inpfdir/$oned_vel .
cp ../../$inpfdir/$data_file .

cp ../../src/prevel3d .
cp ../../src/initinpv .
cp ../../src/migraloc .


# Generate synthetic velocity model from 1D velocity model.
echo "$node_structure" > prevel3d.inp
echo "syndvel.txt" >> prevel3d.inp
echo "synvel.txt" >> prevel3d.inp
echo "synvel_init.txt" >> prevel3d.inp
echo "0.1" >> prevel3d.inp
echo "$oned_vel" >> prevel3d.inp
echo "$topo_file" >> prevel3d.inp
echo "$topo_minlon $topo_maxlon" >> prevel3d.inp
echo "$topo_minlan $topo_maxlan" >> prevel3d.inp
echo "$topo_xnum $topo_ynum" >> prevel3d.inp
echo "$topo_depth" >> prevel3d.inp
echo "$topo_vair" >> prevel3d.inp
./prevel3d

awk '{print $1,$3}' $oned_vel > temp1d_vel.txt
# Generate initial velocity model from 1D velocity model.
echo "init_vel.txt" > initinpv.inp
echo "init_vel_s.txt" >> initinpv.inp
echo "synvel_init.txt" >> initinpv.inp
echo "temp1d_vel.txt" >> initinpv.inp
./initinpv

# Migrate the locations above surface below the surface.
echo "$data_file" >> migraloc.inp
echo "data_temp.txt" >> migraloc.inp
echo "$sta_num" >> migraloc.inp
echo "$topo_file" >> migraloc.inp
echo "$topo_minlon $topo_maxlon" >> migraloc.inp
echo "$topo_minlan $topo_maxlan" >> migraloc.inp
echo "$topo_xnum $topo_ynum" >> migraloc.inp
./migraloc
mv data_temp.txt $data_file

cp  $data_file init_vel.txt ../outf

# Exit working dir
cd ../../
