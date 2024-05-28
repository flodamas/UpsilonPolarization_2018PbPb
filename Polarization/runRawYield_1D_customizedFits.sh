#!/bin/bash

# Define the ROOT environment variables
#ROOTSYS=/path/to/your/root/installation
#source $ROOTSYS/bin/thisroot.sh

# Function to run ROOT script with arguments
run_root_script() {
    local ptMin=$1
    local ptMax=$2
    local refFrameName=$3
    local nCosThetaBins=$4
    local phiMin=$5
    local phiMax=$6
    local iState=$7

    root -l -q -b "rawYield_1D_customizedFits.C(${ptMin}, ${ptMax}, \"${refFrameName}\", ${nCosThetaBins}, ${phiMin}, ${phiMax}, ${iState})"

    #local args=("$@")
    #root -l "rawYield_1D_customizedFits.C(${args[@]})"
}

# Define an array of argument sets
argument_sets=(
#    "0 2 HX 10 -180 180 gUpsilonState"
#    "0 2 HX 8 -180 180 gUpsilonState"
#    "0 2 HX 6 -180 180 gUpsilonState"
#    "2 4 HX 10 -180 180 gUpsilonState"
#    "2 4 HX 8 -180 180 gUpsilonState"
#    "2 4 HX 6 -180 180 gUpsilonState"
#    "4 6 HX 8 -180 180 gUpsilonState"
#    "4 6 HX 6 -180 180 gUpsilonState"
#    "6 8 HX 9 -180 180 gUpsilonState"
#    "6 8 HX 7 -180 180 gUpsilonState"
#    "6 8 HX 5 -180 180 gUpsilonState"
#    "8 12 HX 7 -180 180 gUpsilonState"
#    "8 12 HX 5 -180 180 gUpsilonState"
     "12 16 HX 8 -180 180 gUpsilonState"
     "12 16 HX 6 -180 180 gUpsilonState"
     "16 30 HX 8 -180 180 gUpsilonState"
     "16 30 HX 6 -180 180 gUpsilonState"
)

# Loop through each argument set and run the ROOT script
for args in "${argument_sets[@]}"; do
    run_root_script $args
done

# Compile and run the ROOT script
#root -l "rawYield_1D_customizedFits.C(0 2 \"HX\" 10 -180 180 gUpsilonState)"
