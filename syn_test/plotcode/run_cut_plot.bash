#!/bin/bash

# cut slice n19.4
echo "out_svd_n19.4.txt" > cutslice_reg2.inp
echo "out_svd_n19.4_sper.txt" >> cutslice_reg2.inp
echo "out_svd_n19.4_vper.txt" >> cutslice_reg2.inp
echo "out_svd_ds_ctrl.txt" >> cutslice_reg2.inp
echo "vel.txt" >> cutslice_reg2.inp
echo "vz.1D_inverted" >> cutslice_reg2.inp
echo "204.50 205.20" >> cutslice_reg2.inp
echo "19.40 19.40" >> cutslice_reg2.inp
echo "0.00 25.00" >> cutslice_reg2.inp
echo "0.005 0.005 0.5" >> cutslice_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cutslice_reg2.inp
echo "203.80 205.383" >> cutslice_reg2.inp
echo "18.80 20.383" >> cutslice_reg2.inp
echo "96 96" >> cutslice_reg2.inp
echo "9.673" >> cutslice_reg2.inp
echo "1.5" >> cutslice_reg2.inp

./cutslice_reg2


# cut model resolution slice n19.4
echo "mdres_n19.4.txt" > cut_ressl_reg2.inp
echo "mdres.txt" >> cut_ressl_reg2.inp
echo "vel.txt" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.40 19.40" >> cut_ressl_reg2.inp
echo "0.00 25.00" >> cut_ressl_reg2.inp
echo "0.005 0.005 0.5" >> cut_ressl_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cut_ressl_reg2.inp
echo "203.80 205.383" >> cut_ressl_reg2.inp
echo "18.80 20.383" >> cut_ressl_reg2.inp
echo "96 96" >> cut_ressl_reg2.inp
echo "9.673" >> cut_ressl_reg2.inp
echo "1.5" >> cut_ressl_reg2.inp

./cut_ressl_reg2



# cut dws slice n19.4
echo "dws_n19.4.txt" > cut_ressl_reg2.inp
echo "dws.txt" >> cut_ressl_reg2.inp
echo "vel.txt" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.40 19.40" >> cut_ressl_reg2.inp
echo "0.00 25.00" >> cut_ressl_reg2.inp
echo "0.005 0.005 0.5" >> cut_ressl_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cut_ressl_reg2.inp
echo "203.80 205.383" >> cut_ressl_reg2.inp
echo "18.80 20.383" >> cut_ressl_reg2.inp
echo "96 96" >> cut_ressl_reg2.inp
echo "9.673" >> cut_ressl_reg2.inp
echo "1.5" >> cut_ressl_reg2.inp

./cut_ressl_reg2



# cut rct slice n19.4
echo "rct_n19.4.txt" > cut_ressl_reg2.inp
echo "rct.txt" >> cut_ressl_reg2.inp
echo "vel.txt" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.40 19.40" >> cut_ressl_reg2.inp
echo "0.00 25.00" >> cut_ressl_reg2.inp
echo "0.005 0.005 0.5" >> cut_ressl_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cut_ressl_reg2.inp
echo "203.80 205.383" >> cut_ressl_reg2.inp
echo "18.80 20.383" >> cut_ressl_reg2.inp
echo "96 96" >> cut_ressl_reg2.inp
echo "9.673" >> cut_ressl_reg2.inp
echo "1.5" >> cut_ressl_reg2.inp

./cut_ressl_reg2


# Plot velocity and pertubation of velocity
bash plot.bash
bash plot_per.bash
