#!/bin/bash

# Compile codes
cd src
rm *.o *.mod
make
cd ..

# Compile plot codes
cd src_plot
rm *.o *.mod
make
cd ..
