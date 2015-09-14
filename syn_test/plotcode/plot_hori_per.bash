#!/bin/bash


xitv=0.005
yitv=0.005
echo $xitv $zitv

gmtset FRAME_PEN 1p FRAME_WIDTH 0.12c LABEL_FONT_SIZE 8p LABEL_OFFSET 0.03c ANNOT_FONT_SIZE_SECONDARY 8p ANNOT_FONT_SIZE_PRIMARY 8p HEADER_FONT_SIZE 8p ANNOT_FONT_PRIMARY 4 ANNOT_FONT_SECONDARY 4 LABEL_FONT 4 HEADER_FONT 4 PAGE_ORIENTATION portrait GRID_PEN_PRIMARY 0.5p COLOR_MODEL cmyk PS_COLOR cmyk


# Part 1 1 (x,y)

velfile4=out_svd_z16_vper.txt

awk '{print $1,$2,$4}' "$velfile4" | xyz2grd -R204.5/205.2/19.2/19.7 -G$velfile4.grd -Ddegree/degree/s -I"$xitv"/"$zitv"

makecpt -T-0.05/0.05/0.0005 -I -Cno_green > $velfile4.cpt

grdimage $velfile4.grd -JM10c -C$velfile4.cpt -R204.5/205.2/19.2/19.7 -B0.2/0.2WseN -K -Y5c> $velfile4.ps
psscale -D5c/-1c/6c/0.5ch -O -K -C$velfile4.cpt -B0.02 >> $velfile4.ps


