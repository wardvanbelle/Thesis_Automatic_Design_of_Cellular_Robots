#!/bin/bash

# move to app directory
cd /home/app/

# initialize all needed packages for julia
julia /home/app/initialize.jl 

# run main biobot script
julia /home/app/Biobot_V1.jl

# exit and return ouptput of last script
exit