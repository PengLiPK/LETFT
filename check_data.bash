#!/bin/bash

# Check if the events number are correct or not
pdata1=t_p.txt
sdata1=t_s.txt

pdata2=p_data_initck.txt
sdata2=s_data_initck.txt

awk '{if(NF>5)print $5,$6,$7,$8}' $pdata1 | sort -u > temp1.txt
awk '{if(NF>5)print $5,$6,$7,$8}' $sdata1 | sort -u > temp2.txt

cat temp1.txt temp2.txt | sort -u | sort -k 1n > temp3.txt


awk '{if(NF>5)print $5,$6,$7,$8}' $pdata2 | sort -u > temp4.txt
awk '{if(NF>5)print $5,$6,$7,$8}' $sdata2 | sort -u > temp5.txt

cat temp4.txt temp5.txt | sort -u | sort -k 1n > temp6.txt


paste temp3.txt temp6.txt | awk '{if($1!=$5)print $0}' > temp7.txt
line=$(wc -l temp7.txt | awk '{print $1}')

echo "Numbers of error lines:" $line

if [ "$line" -eq 0 ]
then
	echo "All events are correct!"
	rm temp[1-7].txt
else
	echo $line "events are not correct!"
	echo "Check temp3.txt, temp6.txt and temp7.txt!"
fi
