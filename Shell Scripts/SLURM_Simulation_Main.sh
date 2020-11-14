#!/bin/bash

# Extract Command Line Arguments
noZeroInflate="$1"
RandomEffect="$2"
metadataType="$3"
nPerSubject="$4"
nMetadata="$5"
readDepth="$6"

# Default Values
spikeMicrobes=0.1
if [[ "$metadataType" == *MV* ]]; then
spikeMetadata=0.2
fi
if [[ "$metadataType" == *UV* ]]; then
spikeMetadata=1
fi

# Default Set of Methods
if [[ "$metadataType" == *UVA* ]]; then
if [[ "$RandomEffect" == *F* ]]; then
methods=(limmaVOOM metagenomeSeq metagenomeSeq2 Spearman Spearman.TSS Spearman.CSS Spearman.TMM Spearman.RLE Spearman.CLR LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE ZINB ZINB.CSS ZINB.TMM ZINB.RLE negbin negbin.CSS negbin.TMM negbin.RLE DESeq2 edgeR ZIB)
fi
fi

if [[ "$metadataType" == *UVA* ]]; then
if [[ "$RandomEffect" == *T* ]]; then
methods=(limmaVOOM metagenomeSeq metagenomeSeq2 Spearman Spearman.TSS Spearman.CSS Spearman.TMM Spearman.RLE Spearman.CLR LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE negbin negbin.CSS negbin.TMM negbin.RLE ZIB)
fi
fi

if [[ "$metadataType" == *UVB* ]]; then
if [[ "$RandomEffect" == *F* ]]; then
methods=(limmaVOOM metagenomeSeq metagenomeSeq2 Wilcoxon Wilcoxon.TSS Wilcoxon.CSS Wilcoxon.TMM Wilcoxon.RLE Wilcoxon.CLR LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE ZINB ZINB.CSS ZINB.TMM ZINB.RLE negbin negbin.CSS negbin.TMM negbin.RLE ANCOM DESeq2 edgeR ZIB)
fi
fi

if [[ "$metadataType" == *UVB* ]]; then
if [[ "$RandomEffect" == *T* ]]; then
methods=(limmaVOOM metagenomeSeq metagenomeSeq2 Wilcoxon Wilcoxon.TSS Wilcoxon.CSS Wilcoxon.TMM Wilcoxon.RLE Wilcoxon.CLR LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE negbin negbin.CSS negbin.TMM negbin.RLE ANCOM ZIB)
fi
fi

if [[ "$metadataType" == *MV* ]]; then
if [[ "$RandomEffect" == *F* ]]; then
methods=(limmaVOOM metagenomeSeq LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE ZINB ZINB.CSS ZINB.TMM ZINB.RLE negbin negbin.CSS negbin.TMM negbin.RLE ZIB)
fi
fi

if [[ "$metadataType" == *MV* ]]; then
if [[ "$RandomEffect" == *T* ]]; then
methods=(limmaVOOM metagenomeSeq LM LM.LOG LM.TSS LM.TSS.AST LM.TSS.LOG LM.TSS.LOGIT LM.TMM LM.TMM.LOG LM.CSS LM.CSS.LOG LM.RLE LM.RLE.LOG LM.CLR LM2 LM2.LOG LM2.TSS LM2.TSS.AST LM2.TSS.LOG LM2.TSS.LOGIT LM2.TMM LM2.TMM.LOG LM2.CSS LM2.CSS.LOG LM2.RLE LM2.RLE.LOG LM2.CLR limma limma.LOG limma.TSS limma.TSS.AST limma.TSS.LOG limma.TSS.LOGIT limma.TMM limma.TMM.LOG limma.CSS limma.CSS.LOG limma.RLE limma.RLE.LOG limma.CLR limma2 limma2.LOG limma2.TSS limma2.TSS.AST limma2.TSS.LOG limma2.TSS.LOGIT limma2.TMM limma2.TMM.LOG limma2.CSS limma2.CSS.LOG limma2.RLE limma2.RLE.LOG limma2.CLR CPLM CPLM.TSS CPLM.CSS CPLM.TMM CPLM.RLE negbin negbin.CSS negbin.TMM negbin.RLE ZIB)
fi
fi


# Main LOOP (60 Core Scenarios)
for methodName in "${methods[@]}"; do
for noZeroInflate in ${noZeroInflate}; do
for RandomEffect in ${RandomEffect}; do
for metadataType in ${metadataType}; do
for nSubjects in 10 20 50 100 200; do
for nPerSubject in ${nPerSubject}; do
for nMicrobes in 100 200 500; do
for spikeMicrobes in ${spikeMicrobes}; do
for nMetadata in ${nMetadata}; do
for spikeMetadata in ${spikeMetadata}; do
for effectSize in 1 2 5 10; do
for readDepth in ${readDepth}; do
#
echo "${methodName} ${noZeroInflate} ${RandomEffect} ${metadataType} ${nSubjects} ${nPerSubject} ${nMicrobes} ${spikeMicrobes} ${nMetadata} ${spikeMetadata} ${effectSize} ${readDepth}"
export methodName noZeroInflate RandomEffect metadataType nSubjects nPerSubject nMicrobes spikeMicrobes nMetadata spikeMetadata effectSize readDepth
#
sbatch Simulation_Main.sh
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done
done
done
done
done
done
done
done
