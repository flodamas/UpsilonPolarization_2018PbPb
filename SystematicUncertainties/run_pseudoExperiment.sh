#!/bin/bash
# run_analysis.sh
# Arguments: $1 = job index

cd /afs/cern.ch/user/s/sayoo/private/Polarization/UpsilonPolarization_2018PbPb/SystematicUncertainties

JOB_INDEX=$1
OUTPUT_FILE=output_${JOB_INDEX}.root

# Run ROOT script with output
root -b -q "pseudoExperiment.C(\"$OUTPUT_FILE\", ${JOB_INDEX})"