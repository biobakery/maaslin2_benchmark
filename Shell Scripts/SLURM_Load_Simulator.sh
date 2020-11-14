#!/bin/bash

# Extract Command Line Arguments
readDepth="$1"


# Main Loop (480 Core Scenarios)
for noZeroInflate in FALSE; do
for RandomEffect in FALSE TRUE; do
for metadataType in UVA UVB MVA MVB; do
for nSubjects in 10 20 50 100 200; do
for nMicrobes in 100 200 500; do
for spikeMicrobes in 0.1; do
for effectSize in 1 2 5 10; do
for readDepth in ${readDepth}; do
if [[ "$RandomEffect" == *F* ]]; then
if [[ "$metadataType" == *UV* ]]; then
for nPerSubject in 1; do
for nMetadata in 1; do
for spikeMetadata in 1; do
#
echo "${noZeroInflate} ${RandomEffect} ${metadataType} ${nSubjects} ${nPerSubject} ${nMicrobes} ${spikeMicrobes} ${nMetadata} ${spikeMetadata} ${effectSize} ${readDepth}"
export noZeroInflate RandomEffect metadataType nSubjects nPerSubject nMicrobes spikeMicrobes nMetadata spikeMetadata effectSize readDepth
#
sbatch Load_Simulator.sh
#
sleep 1 # pause to be kind to the scheduler
done
done
done
fi
if [[ "$metadataType" == *MV* ]]; then
for nPerSubject in 1; do
for nMetadata in 5; do
for spikeMetadata in 0.2; do
#
echo "${noZeroInflate} ${RandomEffect} ${metadataType} ${nSubjects} ${nPerSubject} ${nMicrobes} ${spikeMicrobes} ${nMetadata} ${spikeMetadata} ${effectSize} ${readDepth}"
export noZeroInflate RandomEffect metadataType nSubjects nPerSubject nMicrobes spikeMicrobes nMetadata spikeMetadata effectSize readDepth
#
sbatch Load_Simulator.sh
#
sleep 1 # pause to be kind to the scheduler
done
done
done
fi
fi
if [[ "$RandomEffect" == *T* ]]; then
if [[ "$metadataType" == *UV* ]]; then
for nPerSubject in 5; do
for nMetadata in 1; do
for spikeMetadata in 1; do
#
echo "${noZeroInflate} ${RandomEffect} ${metadataType} ${nSubjects} ${nPerSubject} ${nMicrobes} ${spikeMicrobes} ${nMetadata} ${spikeMetadata} ${effectSize} ${readDepth}"
export noZeroInflate RandomEffect metadataType nSubjects nPerSubject nMicrobes spikeMicrobes nMetadata spikeMetadata effectSize readDepth
#
sbatch Load_Simulator.sh
#
sleep 1 # pause to be kind to the scheduler
done
done
done
fi
if [[ "$metadataType" == *MV* ]]; then
for nPerSubject in 5; do
for nMetadata in 5; do
for spikeMetadata in 0.2; do
#
echo "${noZeroInflate} ${RandomEffect} ${metadataType} ${nSubjects} ${nPerSubject} ${nMicrobes} ${spikeMicrobes} ${nMetadata} ${spikeMetadata} ${effectSize} ${readDepth}"
export noZeroInflate RandomEffect metadataType nSubjects nPerSubject nMicrobes spikeMicrobes nMetadata spikeMetadata effectSize readDepth
#
sbatch Load_Simulator.sh
#
sleep 1 # pause to be kind to the scheduler
done
done
done
fi
fi
done
done
done
done
done
done
done
done
