#!/bin/bash


inpf=$1


paste $1 syndata.txt | awk '{if(NF==18)print $5,$6-$15,$7-$16,$8-$17,$9-0.1}' | sort -u > tempdiffloc

num=$(wc -l tempdiffloc | awk '{print $1}')
aadx=$(awk '{a+=sqrt($2*$2)}END{print a/"'$num'"}' tempdiffloc)
aady=$(awk '{a+=sqrt($3*$3)}END{print a/"'$num'"}' tempdiffloc)
aadz=$(awk '{a+=sqrt($4*$4)}END{print a/"'$num'"}' tempdiffloc)
aadt=$(awk '{a+=sqrt($5*$5)}END{print a/"'$num'"}' tempdiffloc)


echo $aadx $aady $aadz $aadt "Average absolute deviations"
echo $aadx $aady $aadz $aadt $inpf "AAD" >> lcresi.txt
