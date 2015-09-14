#!/bin/bash


for((nid=1;nid<=25465;nid++))
do
#	num=$(awk '{if($2=='"$nid"')print $0}' locfdasc.txt | wc -l | awk '{print $1}')
#	echo $num $nid >> statfdloc.txt
	grep " $nid " fd.txt | wc -l | awk '{print $1}' >> temp
done
