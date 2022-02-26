#!/bin/bash

# move to app directory
cd /home/app/

# initialize all needed packages for julia
julia /home/app/initialize.jl 

# install voxcraft and cmake
rm voxcraft-sim -rf; git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;
sudo apt-get update -y; sudo apt-get update; sudo apt-get install -y cmake libboost-all-dev
cd voxcraft-sim; mkdir build; cd build; cmake ..; make -j 10

# run main biobot script
julia /home/app/Biobot_V1.jl

# exit and return ouptput of last script
exit