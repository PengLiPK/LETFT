#!/bin/bash


tmpresl=$1
resln=+$tmpresl
xyint=0.5/0.5
xyrfint=0.1/0.1

gmt gmtset MAP_TICK_LENGTH_PRIMARY 2.5p/1.2p MAP_ANNOT_OFFSET_PRIMARY 2p MAP_FRAME_PEN 1p MAP_FRAME_WIDTH 0.12c FONT_LABEL 8p MAP_LABEL_OFFSET 0.2c FONT_ANNOT_SECONDARY 8p FONT_ANNOT_PRIMARY 8p FONT_TITLE 8p PS_PAGE_ORIENTATION portrait MAP_GRID_PEN_PRIMARY 0.5p COLOR_MODEL cmyk


vf11=slice_veltp_1.txt
resf1=slice_mdres_1.txt
elel1=elev_line1.txt
elep1=elev_poly1.txt
ndprj1=ndprj_1.txt

outf=veltp_line
#xmax=$(tail -1 $vf11 | awk '{print $1}')

gmt makecpt -T2.7/8.2/0.01 -D -I -Cwysiwyg > $outf.cpt


#################################################################################################
# Part 1 1 (x,y)
xmax=$(tail -1 $elel1 |awk '{printf "%d",$1}')

awk '{print $1,-4.079+$2,$3}' "$resf1" | gmt surface -R0/$xmax/-5/25 -G$resf1.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $resf1.nm.grd -G$resf1.rf.grd -I$xyrfint
awk '{print $1,-4.079+$2}' $elep1 | gmt grdmask -G$elep1.grd -N1/1/NaN -R -I$xyrfint
gmt grdmath $resf1.rf.grd $elep1.grd OR = $resf1.grd


# File for mask data 
gmt grd2xyz $resf1.rf.grd > temp1.txt
num=$(wc -l temp1.txt | awk '{print $1}')

awk '{if(NF==3)print $1,-4.079+$2,$3}' "$vf11" | gmt surface -R0/$xmax/-5/25 -G$vf11.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $vf11.nm.grd -G$vf11.rf.grd -I$xyrfint
gmt grd2xyz $vf11.rf.grd > temp2.txt

./valmask temp1.txt temp2.txt $vf11.mk.txt $tmpresl

gmt xyz2grd $vf11.mk.txt -G$vf11.mk.grd  -R0/$xmax/-5/25 -I$xyrfint
gmt grdmath $vf11.mk.grd $elep1.grd OR = $vf11.grd
gmt grdimage $vf11.grd -JX9.6c/-4c -C$outf.cpt -R0/$xmax/-4.079/21 -B5:Distance\ \(km\):/5:Depth\ \(km\):WSen -K > $outf.ps
awk '{print $1,-4.079+$2}' $elel1 | gmt psxy -J -R -W2p -O -K >>$outf.ps
#awk '{print $1,-4.079+$2}' nodelys.txt | gmt psxy -J -R -W1p,- -O -K >> $outf.ps

# Plot lines
for x in $(awk '{if(NF==1) print $1}' $vf11)
do
	echo $x -4.079 > templine
	echo $x 21.0 >> templine
	gmt psxy templine -J -R -W1p,- -O -K >> $outf.ps
done

awk '{print $1,-4.079+$2}' $ndprj1 | gmt psxy -J -R -S+0.15 -Gwhite -W1.2p,white -K -O >> $outf.ps

##############################################################################################

gmt psscale -D-3.5c/-1.2c/10c/0.3ch -O -C$outf.cpt -B1:Velocity\ \(km\/s\): >> $outf.ps


gmt ps2raster -A -Te $outf.ps
