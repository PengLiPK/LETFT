#!/bin/bash


resf=$1


num=$(wc -l $resf | awk '{print $1}')
mres=$(awk '{a+=$1}END{print a/"'$num'"}' $resf)
stdres=$(awk '{a+=($1-"'$mres'")*($1-"'$mres'")}END{print sqrt(a/"'$num'")}' $resf)
advres=$(awk '{a+=sqrt($1*$1)}END{print a/"'$num'"}' $resf)

echo $mres $stdres $advres $resf

echo $mres $stdres $advres $resf >> statres.txt
