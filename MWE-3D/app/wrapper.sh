#!/bin/bash

# move to app directory
cd /home/app/

# initialize all needed packages for julia
julia /home/app/initialize.jl 

# install voxcraft and cmake
rm voxcraft-sim -rf; git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;
apt-get update; apt-get install -y cmake

# run main biobot script
julia /home/app/Biobot_V1.jl

# exit and return ouptput of last script
exit