#!/bin/bash

gmt pshistogram -Bxa0.02f0.1+l"Noise" -Bya2f1+l"Frequency"+u" %" \
	-BWSne+t"Histograms"+glightblue noiseGS.txt -R-0.1/0.1/0/9 \
	-JX10c/10c -Gorange -L1p -Z1 -W0.001 -P > noiseGS.ps
