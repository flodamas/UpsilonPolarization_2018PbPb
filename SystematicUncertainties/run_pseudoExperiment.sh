#!/bin/bash
# run_analysis.sh
# Arguments: $1 = job index

cd /afs/cern.ch/user/s/sayoo/private/Polarization/UpsilonPolarization_2018PbPb/SystematicUncertainties

JOB_INDEX=$1

# Number of bins for cosTheta and phi
cosThetaBins=5
phiBins=3

# Define ranges for cosThetaPhi and phi bins
cosThetaMin=-0.7
cosThetaMax=0.7
cosThetaStep=$(echo "($cosThetaMax - $cosThetaMin) / $cosThetaBins" | bc -l)

# Define ranges for phi bins
phiMin=0
phiMax=180
phiStep=$(echo "($phiMax - $phiMin) / $phiBins" | bc -l)

# Loop over cosTheta and phi bins
for i in $(seq 0 $((cosThetaBins - 1))); do
    cosThetaLow=$(echo "$cosThetaMin + $i * $cosThetaStep" | bc -l)
    cosThetaHigh=$(echo "$cosThetaLow + $cosThetaStep" | bc -l)

    for j in $(seq 0 $((phiBins - 1))); do
        phiLow=$(echo "$phiMin + $j * $phiStep" | bc -l)
        phiHigh=$(echo "$phiLow + $phiStep" | bc -l)

        # Format cosTheta and phi values to two decimal places
        cosThetaLowFormatted=$(printf "%.2f" $cosThetaLow)
        cosThetaHighFormatted=$(printf "%.2f" $cosThetaHigh)
        phiLowFormatted=$(printf "%.2f" $phiLow)
        phiHighFormatted=$(printf "%.2f" $phiHigh)

        # Define the output file with formatted values
        OUTPUT_FILE="output_${JOB_INDEX}_cosTheta_${cosThetaLowFormatted}_${cosThetaHighFormatted}_phi_${phiLowFormatted}_${phiHighFormatted}.root"

        # Run ROOT script with output
        root -b -q "pseudoExperiment_condor.C(\"$OUTPUT_FILE\", ${JOB_INDEX}, 0, 2, kFALSE, $cosThetaLow, $cosThetaHigh, $phiLow, $phiHigh)"
    done
done