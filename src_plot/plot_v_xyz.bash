#!/bin/bash


tmpresl=$1
resln=+$tmpresl
xyint=0.5/0.5
xyrfint=0.1/0.1

gmt gmtset MAP_TICK_LENGTH_PRIMARY 2.5p/1.2p MAP_ANNOT_OFFSET_PRIMARY 2p MAP_FRAME_PEN 1p MAP_FRAME_WIDTH 0.12c FONT_LABEL 8p MAP_LABEL_OFFSET 5p FONT_ANNOT_SECONDARY 8p FONT_ANNOT_PRIMARY 8p FONT_TITLE 8p PS_PAGE_ORIENTATION portrait MAP_GRID_PEN_PRIMARY 0.5p COLOR_MODEL cmyk MAP_ANNOT_OFFSET_PRIMARY 2p MAP_TICK_LENGTH_PRIMARY 2p/1p PS_COLOR_MODEL cmyk


vf11=slice_truevel.txt
vf12=slice_veltp.txt
vf13=slice_veldifftp.txt
resf=slice_mdres.txt
initvf=vz.1D_inverted

outf=syn_comp_xyz
#xmax=$(tail -1 $vf11 | awk '{print $1}')
xmax=70

awk '{print $1,-4.079+$3,$4}' "$resf" | gmt surface -R0/$xmax/-5/25 -G$resf.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $resf.nm.grd -G$resf.rf.grd -I$xyrfint
awk '{print $1,-4.079+$3}' elev_poly.txt | gmt grdmask -Gelev_poly.grd -N1/1/NaN -R -I$xyrfint
gmt grdmath $resf.rf.grd elev_poly.grd OR = $resf.grd

gmt makecpt -T-12/12/0.01 -D -I -Cpolar > $outf.cpt


# Part 1 1 (1-D)
initvnum=$(head -1 $initvf | awk '{print $1}')
echo $initvnum $initvf
tail -$initvnum $initvf | awk '{print $3, -4.079+$1}' | gmt psxy -JX2.6c/-2.6c -R0/9/-4.079/21 -W1p -B2:Velocity\ \(km\/s\):/5:Depth\ \(km\):WseN -K -X2.5c -Y10c > $outf.ps


# Part 2 1 (x,y)
awk '{print $1,-4.079+$3,100*$4}' "$vf11" | gmt surface -R0/$xmax/-5/25 -G$vf11.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $vf11.nm.grd -G$vf11.rf.grd -I$xyrfint
gmt grdmath $vf11.rf.grd elev_poly.grd OR = $vf11.grd

gmt grdimage $vf11.grd -JX5c/-2.6c -C$outf.cpt -R0/$xmax/-4.079/21 -B10/5:Depth\ \(km\):WseN -K -O -X-1c -Y-3.5c >> $outf.ps
awk '{print $1,-4.079+$3}' elev_line.txt | gmt psxy -J -R -W2p -O -K >>$outf.ps


# Part 1 2 (x,y)
awk '{print $1,-4.079+$3,100*$4}' "$vf12" | gmt surface -R0/$xmax/-5/25 -G$vf12.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $vf12.nm.grd -G$vf12.rf.grd -I$xyrfint
gmt grdmath $vf12.rf.grd elev_poly.grd OR = $vf12.grd

gmt grdimage $vf12.grd -JX5c/-2.6c -C$outf.cpt -R0/$xmax/-4.079/21 -B10:Distance\ \(km\):/5:Depth\ \(km\):WseN -O -K -X6c -Y3.5c>> $outf.ps
gmt grdcontour $resf.grd -J -R -C$resln -W1p,black -S4 -O -K >> $outf.ps
awk '{print $1,-4.079+$3}' elev_line.txt | gmt psxy -J -R -W2p -O -K >>$outf.ps


# File for mask data 
gmt grd2xyz $resf.rf.grd > temp1.txt
num=$(wc -l temp1.txt | awk '{print $1}')

# Part 2 2 (x,y)
awk '{print $1,-4.079+$3,100*$4}' "$vf13" | gmt surface -R0/$xmax/-5/25 -G$vf13.nm.grd -I$xyint -T0.25 -C0.1
gmt grdsample $vf13.nm.grd -G$vf13.rf.grd -I$xyrfint
gmt grd2xyz $vf13.rf.grd > temp2.txt

./valmask temp1.txt temp2.txt $vf13.mk.txt $tmpresl

gmt xyz2grd $vf13.mk.txt -G$vf13.mk.grd  -R0/$xmax/-5/25 -I$xyrfint
gmt grdmath $vf13.mk.grd elev_poly.grd OR = $vf13.grd
gmt grdimage $vf13.grd -JX5c/-2.6c -C$outf.cpt -R0/$xmax/-4.079/21 -B10/5WseN -O -K -Y-3.5c>> $outf.ps
awk '{print $1,-4.079+$3}' elev_line.txt | gmt psxy -J -R -W2p -O -K >>$outf.ps


gmt psscale -D-0.9c/-0.5c/6c/0.3ch -O -K -C$outf.cpt -B4:Velocity\ perturbation\ \(\%\): >> $outf.ps


gmt pstext -R0/20/0/20 -Jx1c -N -O -X-15c -Y-5c << EOF >> $outf.ps
8.5 11.5 8 0.0 0 CB a)
14.5.5 11.5 8 0.0 0 CB c)
8.5 7.9 8 0.0 0 CB b)
14.5 7.9 8 0.0 0 CB d)
EOF

gmt ps2raster -A -Te $outf.ps
