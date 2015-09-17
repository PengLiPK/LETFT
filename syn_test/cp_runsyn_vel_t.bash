#!/bin/bash

wkdir = work
inpfdir = input_file

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
cp ../../$input_file/vstrct_tp3km.txt .
cp ../../$input_file/hawaii_tp4tomo2.xyz .
cp ../../$input_file/vz.1D_inverted .
