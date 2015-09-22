#!/bin/bash

inpf=$1


# cut slice n19.4
echo "out_svd_n19.4v.txt" > cutv_reg2.inp
echo "$inpf" >> cutv_reg2.inp
echo "204.50 205.20" >> cutv_reg2.inp
echo "19.41 19.41" >> cutv_reg2.inp
echo "0.00 25.00" >> cutv_reg2.inp
echo "0.005 0.005 0.5" >> cutv_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cutv_reg2.inp
echo "203.80 205.383" >> cutv_reg2.inp
echo "18.80 20.383" >> cutv_reg2.inp
echo "96 96" >> cutv_reg2.inp
echo "9.673" >> cutv_reg2.inp
echo "1.5" >> cutv_reg2.inp

./cutv_reg2



# cut slice n19.4
echo "out_svd_n19.4vdiff.txt" > cutvdiff_reg2.inp
echo "mddiff.txt" >> cutvdiff_reg2.inp
echo "mdres.txt 0.3" >> cutvdiff_reg2.inp
echo "204.50 205.20" >> cutvdiff_reg2.inp
echo "19.41 19.41" >> cutvdiff_reg2.inp
echo "0.00 25.00" >> cutvdiff_reg2.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg2.inp
echo "203.80 205.383" >> cutvdiff_reg2.inp
echo "18.80 20.383" >> cutvdiff_reg2.inp
echo "96 96" >> cutvdiff_reg2.inp
echo "9.673" >> cutvdiff_reg2.inp
echo "1.5" >> cutvdiff_reg2.inp

./cutvdiff_reg2


# cut slice n19.4
echo "out_svd_n19.4vdiffper.txt" > cutvdiff_reg2.inp
echo "mddiffper.txt" >> cutvdiff_reg2.inp
echo "mdres.txt 0.3" >> cutvdiff_reg2.inp
echo "204.50 205.20" >> cutvdiff_reg2.inp
echo "19.41 19.41" >> cutvdiff_reg2.inp
echo "0.00 25.00" >> cutvdiff_reg2.inp
echo "0.005 0.005 0.5" >> cutvdiff_reg2.inp
echo "hawaii_tp4tomo2.xyz" >> cutvdiff_reg2.inp
echo "203.80 205.383" >> cutvdiff_reg2.inp
echo "18.80 20.383" >> cutvdiff_reg2.inp
echo "96 96" >> cutvdiff_reg2.inp
echo "9.673" >> cutvdiff_reg2.inp
echo "1.5" >> cutvdiff_reg2.inp

./cutvdiff_reg2


# cut model resolution slice n19.4
echo "mdres_n19.4.txt" > cut_ressl_reg2.inp
echo "mdres.txt" >> cut_ressl_reg2.inp
echo "$inpf" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.41 19.41" >> cut_ressl_reg2.inp
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
echo "$inpf" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.41 19.41" >> cut_ressl_reg2.inp
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
echo "$inpf" >> cut_ressl_reg2.inp
echo "204.50 205.20" >> cut_ressl_reg2.inp
echo "19.41 19.41" >> cut_ressl_reg2.inp
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
bash plot_v.bash
bash plot_vdiff.bash
bash plot_vdiffper.bash
