#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=work
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

# Parameters of synthetic velocity model construction
#######################################################
node_structure=vstrct_tp3km.txt
oned_vel=vz.1D_inverted
vel_pertubation=0.1 # 0.1 is 10% from 1D velocity model.
#######################################################



# Check working dir is exist or not.
if [ ! -d $wkdir ]
then
	mkdir $wkdir
fi

if [ ! -d $wkdir/syn_vel_t/ ]
then
	mkdir $wkdir/syn_vel_t/
fi

# Enter $wkdir/syn_vel_t
cd $wkdir/syn_vel_t

# Copy input files
cp ../../$inpfdir/$node_structure .
cp ../../$inpfdir/$topo_file .
cp ../../$inpfdir/$oned_vel .

cp ../../src/prevel3d .
cp ../../src/fmm_synt3d_v .

# Generate synthetic velocity model from 1D velocity model.
echo "$node_structure" > prevel3d.inp
echo "syndvel_10per.txt" >> prevel3d.inp
echo "synvel_10per.txt" >> prevel3d.inp
echo "synvel_init.txt" >> prevel3d.inp
echo "$vel_pertubation" >> prevel3d.inp
echo "$oned_vel" >> prevel3d.inp
echo "$topo_file" >> prevel3d.inp
echo "$topo_minlon $topo_maxlon" >> prevel3d.inp
echo "$topo_minlan $topo_maxlan" >> prevel3d.inp
echo "$topo_xnum $topo_ynum" >> prevel3d.inp
echo "$topo_depth" >> prevel3d.inp
echo "$topo_vair" >> prevel3d.inp
./prevel3d

# Generate synthetic travel time data from synthetic velocity model.
echo "$data_file" > fmm_synt3d_v.inp
echo "synvel_10per.txt" >> fmm_synt3d_v.inp
echo "$sta_num" >> fmm_synt3d_v.inp
echo "$fmm_dx $fmm_dy $fmm_dz" >> fmm_synt3d_v.inp
echo "$vel_pertubation" >> fmm_synt3d_v.inp
echo "$oned_vel" >> fmm_synt3d_v.inp
echo "$topo_file" >> fmm_synt3d_v.inp
echo "$topo_minlon $topo_maxlon" >> fmm_synt3d_v.inp
echo "$topo_minlan $topo_maxlan" >> fmm_synt3d_v.inp
echo "$topo_xnum $topo_ynum" >> fmm_synt3d_v.inp
echo "$topo_depth" >> fmm_synt3d_v.inp
echo "$topo_vair" >> fmm_synt3d_v.inp
./fmm_synt3d_v
