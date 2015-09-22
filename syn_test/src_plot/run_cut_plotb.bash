#!/bin/bash

# cut slice n19.4
echo "out_svd_n19.4.txt" > cutslice_reg2b.inp
echo "out_svd_n19.4_sper.txt" >> cutslice_reg2b.inp
echo "out_svd_n19.4_vper.txt" >> cutslice_reg2b.inp
echo "out_svd_ds.txt" >> cutslice_reg2b.inp
echo "vel.txt" >> cutslice_reg2b.inp
echo "vz.1D_inverted" >> cutslice_reg2b.inp
echo "204.50 205.20" >> cutslice_reg2b.inp
echo "19.40 19.40" >> cutslice_reg2b.inp
echo "0.00 55.00" >> cutslice_reg2b.inp
echo "0.005 0.005 0.5" >> cutslice_reg2b.inp
echo "hawaii_tp4tomo2.xyz" >> cutslice_reg2b.inp
echo "203.80 205.383" >> cutslice_reg2b.inp
echo "18.80 20.383" >> cutslice_reg2b.inp
echo "96 96" >> cutslice_reg2b.inp
echo "9.673" >> cutslice_reg2b.inp
echo "1.5" >> cutslice_reg2b.inp


./cutslice_reg2b



# Plot velocity and pertubation of velocity
bash plot.bash
bash plot_perb.bash
