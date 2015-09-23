#!/bin/bash

# Parameters of working dir
#######################################################
wkdir=syn_work_ps
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
topo_vair_p=2.0
topo_vair_s=1.15
#######################################################

# Parameters of synthetic velocity model construction
#######################################################
node_structure=vstrct_tp3km.txt
oned_vel_p=vz.1D_inverted
oned_vel_s=vz.1D_inverted_s
vel_pertubation=0.1 # 0.1 is 10% from 1D velocity model.
#######################################################

# Parameters of synthetic travel times construction
#######################################################
data_file_p=p_data.txt
data_file_s=s_data.txt
sta_num_p=35
sta_num_s=35
minlon=204.50
maxlon=205.20
minlan=19.20
maxlan=19.70
minz=0
maxz=25.0
fmm_dx=0.005
fmm_dy=0.005
fmm_dz=0.5
fmm_rf_nx=5.0
fmm_rf_ny=5.0
fmm_rf_nz=5.0
fmm_rf_x=0.045
fmm_rf_y=0.045
fmm_rf_z=5
#######################################################

# Paramters of constructing input data for synthetic tomograpy.
#######################################################
nsmean=0
nslcstd=1
nsotstd=0.1
nevn=1800
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

if [ ! -d $wkdir/syn_outf/ ]
then
	mkdir $wkdir/syn_outf/
fi

# Enter $wkdir/syn_vel_t
cd $wkdir/syn_vel_t

# Copy input files
cp ../../$inpfdir/$node_structure .
cp ../../$inpfdir/$topo_file .
cp ../../$inpfdir/$oned_vel .
cp ../../$inpfdir/$data_file .

cp ../../src/prevel3d .
cp ../../src/fmm_synt3d_v .
cp ../../src/genGSnoise .
cp ../../src/migraloc .

# Generate synthetic P-velocity model from 1D P-velocity model.
echo "$node_structure" > prevel3d.inp
echo "syndvelp.txt" >> prevel3d.inp
echo "synvelp.txt" >> prevel3d.inp
echo "synvelp_init.txt" >> prevel3d.inp
echo "$vel_pertubation" >> prevel3d.inp
echo "$oned_vel_p" >> prevel3d.inp
echo "$topo_file" >> prevel3d.inp
echo "$topo_minlon $topo_maxlon" >> prevel3d.inp
echo "$topo_minlan $topo_maxlan" >> prevel3d.inp
echo "$topo_xnum $topo_ynum" >> prevel3d.inp
echo "$topo_depth" >> prevel3d.inp
echo "$topo_vair_p" >> prevel3d.inp
./prevel3d

# Generate synthetic S-velocity model from 1D S-velocity model.
echo "$node_structure" > prevel3d.inp
echo "syndvels.txt" >> prevel3d.inp
echo "synvels.txt" >> prevel3d.inp
echo "synvels_init.txt" >> prevel3d.inp
echo "$vel_pertubation" >> prevel3d.inp
echo "$oned_vel_s" >> prevel3d.inp
echo "$topo_file" >> prevel3d.inp
echo "$topo_minlon $topo_maxlon" >> prevel3d.inp
echo "$topo_minlan $topo_maxlan" >> prevel3d.inp
echo "$topo_xnum $topo_ynum" >> prevel3d.inp
echo "$topo_depth" >> prevel3d.inp
echo "$topo_vair_s" >> prevel3d.inp
./prevel3d

# Generate synthetic P-travel time data from synthetic P-velocity model.
echo "$data_file_p" > fmm_synt3d_v.inp
echo "synvelp.txt" >> fmm_synt3d_v.inp
echo "$sta_num_p" >> fmm_synt3d_v.inp
echo "$fmm_dx $fmm_dy $fmm_dz" >> fmm_synt3d_v.inp
echo "$fmm_rf_nx $fmm_rf_ny $fmm_rf_nz" >> fmm_synt3d_v.inp
echo "$fmm_rf_x $fmm_rf_y $fmm_rf_z" >> fmm_synt3d_v.inp
echo "$minlon $maxlon" >> fmm_synt3d_v.inp
echo "$minlan $maxlan" >> fmm_synt3d_v.inp
echo "$minz $maxz" >> fmm_synt3d_v.inp
echo "2" >> fmm_synt3d_v.inp
echo "2" >> fmm_synt3d_v.inp
echo "$topo_file" >> fmm_synt3d_v.inp
echo "$topo_minlon $topo_maxlon" >> fmm_synt3d_v.inp
echo "$topo_minlan $topo_maxlan" >> fmm_synt3d_v.inp
echo "$topo_xnum $topo_ynum" >> fmm_synt3d_v.inp
echo "$topo_depth" >> fmm_synt3d_v.inp
echo "$topo_vair_p" >> fmm_synt3d_v.inp
./fmm_synt3d_v
mv t.txt syndata_p.txt

# Generate synthetic P-travel time data from synthetic P-velocity model.
echo "$data_file_s" > fmm_synt3d_v.inp
echo "synvels.txt" >> fmm_synt3d_v.inp
echo "$sta_num_s" >> fmm_synt3d_v.inp
echo "$fmm_dx $fmm_dy $fmm_dz" >> fmm_synt3d_v.inp
echo "$fmm_rf_nx $fmm_rf_ny $fmm_rf_nz" >> fmm_synt3d_v.inp
echo "$fmm_rf_x $fmm_rf_y $fmm_rf_z" >> fmm_synt3d_v.inp
echo "$minlon $maxlon" >> fmm_synt3d_v.inp
echo "$minlan $maxlan" >> fmm_synt3d_v.inp
echo "$minz $maxz" >> fmm_synt3d_v.inp
echo "2" >> fmm_synt3d_v.inp
echo "2" >> fmm_synt3d_v.inp
echo "$topo_file" >> fmm_synt3d_v.inp
echo "$topo_minlon $topo_maxlon" >> fmm_synt3d_v.inp
echo "$topo_minlan $topo_maxlan" >> fmm_synt3d_v.inp
echo "$topo_xnum $topo_ynum" >> fmm_synt3d_v.inp
echo "$topo_depth" >> fmm_synt3d_v.inp
echo "$topo_vair_s" >> fmm_synt3d_v.inp
./fmm_synt3d_v
mv t.txt syndata_s.txt

# Add Gaussian noise to the locations and origin times of synthetic data.
# The output files will be used as input data for synthetic tomography.
echo "minlon = $minlon" > addgns4loc_config.py
echo "maxlon = $maxlon" >> addgns4loc_config.py
echo "minlan = $minlan" >> addgns4loc_config.py
echo "maxlan = $maxlan" >> addgns4loc_config.py
echo "minz = $minz" >> addgns4loc_config.py
echo "maxz = $maxz" >> addgns4loc_config.py
echo "nsmean = $nsmean" >> addgns4loc_config.py
echo "nslcstd = $nslcstd" >> addgns4loc_config.py
echo "nsotstd = $nsotstd" >> addgns4loc_config.py
echo "nevn = $nevn" >> addgns4loc_config.py
echo "nstap = $sta_num_p" >> addgns4loc_config.py
echo "inpfp = 'syndata_p.txt'" >> addgns4loc_config.py
echo "outfp = 'syndat_p_init_temp.txt'" >> addgns4loc_config.py
echo "nstas = $sta_num_s" >> addgns4loc_config.py
echo "inpfs = 'syndata_s.txt'" >> addgns4loc_config.py
echo "outfs = 'syndata_s_init_temp.txt'" >> addgns4loc_config.py
python addgns4loc.py

# Migrate the locations above surface below the surface for P wave data.
echo "syndata_p_init_temp.txt" >> migraloc.inp
echo "syndata_P_init.txt" >> migraloc.inp
echo "$sta_num_p" >> migraloc.inp
echo "$topo_file" >> migraloc.inp
echo "$topo_minlon $topo_maxlon" >> migraloc.inp
echo "$topo_minlan $topo_maxlan" >> migraloc.inp
echo "$topo_xnum $topo_ynum" >> migraloc.inp
./migraloc

# Migrate the locations above surface below the surface for S wave data.
echo "syndata_s_init_temp.txt" >> migraloc.inp
echo "syndata_s_init.txt" >> migraloc.inp
echo "$sta_num_p" >> migraloc.inp
echo "$topo_file" >> migraloc.inp
echo "$topo_minlon $topo_maxlon" >> migraloc.inp
echo "$topo_minlan $topo_maxlan" >> migraloc.inp
echo "$topo_xnum $topo_ynum" >> migraloc.inp
./migraloc

cp syndata_p.txt syndata_p_init.txt synvelp.txt synvelp_init.txt ../syn_outf
cp syndata_s.txt syndata_s_init.txt synvels.txt synvels_init.txt ../syn_outf

# Exit working dir
cd ../../
