#!/bin/bash



awk '{if($4==19.40000000) print $3,$4,$5,$2}' vfile.txt > n19.4_vper.txt
#awk '{if($2==19.41000000000000) print $0}' synvelper.txt > n19.41_vper.txt
