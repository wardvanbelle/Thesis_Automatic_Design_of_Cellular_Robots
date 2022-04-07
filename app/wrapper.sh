#!/bin/bash

# move to app directory
cd Thesis_Automatic_Design_of_Cellular_Robots/app/

# initialize all needed packages for julia
julia initialize.jl 

# install voxcraft and cmake
apt update;apt install -y cmake libboost-all-dev;
rm voxcraft-sim -rf; git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;
mkdir build; cd build; cmake .. -DBOOST_INCLUDEDIR="/usr/include" -DBOOST_LIBRARYDIR="/usr/lib/x86_64-linux-gnu"; ln -s /usr/include/boost /opt/conda/bin/../x86_64-conda-linux-gnu/sysroot/usr/include/boost; make -j 10;

# run main biobot script
cd ../../
julia Biobot_V1.jl

# move the files of the best result to the project folder
mv ./Biobot_V1/histories/best* /project
mv ./Biobot_V1/xmls/best* /project

# exit and return ouptput of last script
exit