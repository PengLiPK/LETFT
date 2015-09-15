#!/bin/bash

# Compile modules
cd modules
rm *.o rm *.mod
make
cd ..

# Compile synthetic velocity and travel time programs
cd syn_t_vel/
rm *.o rm *.mod
cp ../modules/*.mod ./
cp ../modules/*.o ./
make
cd ..