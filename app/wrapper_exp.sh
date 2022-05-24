#!/bin/bash

# move to app directory
cd Thesis_Automatic_Design_of_Cellular_Robots/app/

# initialize all needed packages for julia en python
julia initialize.jl 
pip install lxml
pip install numpy

# install voxcraft and cmake
apt update;apt install -y cmake libboost-all-dev zip;
rm voxcraft-sim -rf; git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;
mkdir build; cd build; cmake .. -DBOOST_INCLUDEDIR="/usr/include" -DBOOST_LIBRARYDIR="/usr/lib/x86_64-linux-gnu"; ln -s /usr/include/boost /opt/conda/bin/../x86_64-conda-linux-gnu/sysroot/usr/include/boost; make -j 10;

# run main biobot script
cd ../../
julia Biobot_exp.jl $1 $2 $3

# exit and return ouptput of last script
exit