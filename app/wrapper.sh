#!/bin/bash

# move to app directory
cd /home/app/

# initialize all needed packages for julia
julia initialize.jl 

# install voxcraft and cmake
sudo apt update; sudo apt install -y cmake libboost-all-dev;
rm voxcraft-sim -rf; git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;
mkdir build; cd build; cmake .. -DBOOST_INCLUDEDIR="/usr/include" -DBOOST_LIBRARYDIR="/usr/lib/x86_64-linux-gnu"; ln -s /usr/include/boost /opt/conda/bin/../x86_64-conda-linux-gnu/sysroot/usr/include/boost; make -j 10;

# run main biobot script
cd ../../
julia Biobot_V1.jl

# exit and return ouptput of last script
exit