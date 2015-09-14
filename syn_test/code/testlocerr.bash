#!/bin/bash


inp1=$1

echo $inp1

awk '{if(($4*$4) <"'$inp1'" )print $0}' tempdiffloc > tmp1

num=$(wc -l tmp1 | awk '{print $1}')
aadx=$(awk '{a+=sqrt($2*$2)}END{print a/"'$num'"}' tmp1)
aady=$(awk '{a+=sqrt($3*$3)}END{print a/"'$num'"}' tmp1)
aadz=$(awk '{a+=sqrt($4*$4)}END{print a/"'$num'"}' tmp1)
aadt=$(awk '{a+=sqrt($5*$5)}END{print a/"'$num'"}' tmp1)


echo $aadx $aady $aadz $aadt "Average absolute deviations"
