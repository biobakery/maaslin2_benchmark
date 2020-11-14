#!/bin/bash
#
#SBATCH --job-name=Simulation_Main
#SBATCH --ntasks=1 # Ensure that you have only one task running
#SBATCH --cpus-per-task=4 # Number of cores to the above tasks
#SBATCH --partition=serial_requeue # Partition to submit to
#
# Time format = HH:MM:SS, DD-HH:MM:SS
#
#SBATCH --time=50:00:00
#
# Minimum memory required per allocated  CPU  in  MegaBytes.
#
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-type=NONE
#SBATCH --mail-user=hmallick@hsph.harvard.edu

module load centos6/0.0.1-fasrc01
source new-modules.sh
module load gcc/5.2.0-fasrc01 libxml2/2.7.8-fasrc02 R/3.3.1-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER

if [[ "$noZeroInflate" == *F* && "$RandomEffect" == *F* ]]; then
Rscript --no-save --no-restore --verbose Simulation_Main.R --methodName ${methodName} --metadataType ${metadataType} --nSubjects ${nSubjects} --nPerSubject ${nPerSubject} --nMicrobes ${nMicrobes} --spikeMicrobes ${spikeMicrobes} --nMetadata ${nMetadata} --spikeMetadata ${spikeMetadata} --effectSize ${effectSize} --readDepth ${readDepth} --workingDirectory /n/regal/huttenhower_lab/hmallick/Maaslin2 > ${methodName}_ZeroInflate_noRandomEffect_${metadataType}_${nSubjects}_${nPerSubject}_${nMicrobes}_${spikeMicrobes}_${nMetadata}_${spikeMetadata}_${effectSize}_${readDepth}.Rout
fi
if [[ "$noZeroInflate" == *T* && "$RandomEffect" == *F* ]]; then
Rscript --no-save --no-restore --verbose Simulation_Main.R --methodName ${methodName} --metadataType ${metadataType} --nSubjects ${nSubjects} --nPerSubject ${nPerSubject} --nMicrobes ${nMicrobes} --spikeMicrobes ${spikeMicrobes} --nMetadata ${nMetadata} --spikeMetadata ${spikeMetadata} --effectSize ${effectSize} --readDepth ${readDepth} --workingDirectory /n/regal/huttenhower_lab/hmallick/Maaslin2 --noZeroInflate > ${methodName}_noZeroInflate_noRandomEffect_${metadataType}_${nSubjects}_${nPerSubject}_${nMicrobes}_${spikeMicrobes}_${nMetadata}_${spikeMetadata}_${effectSize}_${readDepth}.Rout
fi
if [[ "$noZeroInflate" == *F* && "$RandomEffect" == *T* ]]; then
Rscript --no-save --no-restore --verbose Simulation_Main.R --methodName ${methodName} --metadataType ${metadataType} --nSubjects ${nSubjects} --nPerSubject ${nPerSubject} --nMicrobes ${nMicrobes} --spikeMicrobes ${spikeMicrobes} --nMetadata ${nMetadata} --spikeMetadata ${spikeMetadata} --effectSize ${effectSize} --readDepth ${readDepth} --workingDirectory /n/regal/huttenhower_lab/hmallick/Maaslin2 --RandomEffect > ${methodName}_ZeroInflate_RandomEffect_${metadataType}_${nSubjects}_${nPerSubject}_${nMicrobes}_${spikeMicrobes}_${nMetadata}_${spikeMetadata}_${effectSize}_${readDepth}.Rout
fi
if [[ "$noZeroInflate" == *T* && "$RandomEffect" == *T* ]]; then
Rscript --no-save --no-restore --verbose Simulation_Main.R --methodName ${methodName} --metadataType ${metadataType} --nSubjects ${nSubjects} --nPerSubject ${nPerSubject} --nMicrobes ${nMicrobes} --spikeMicrobes ${spikeMicrobes} --nMetadata ${nMetadata} --spikeMetadata ${spikeMetadata} --effectSize ${effectSize} --readDepth ${readDepth} --workingDirectory /n/regal/huttenhower_lab/hmallick/Maaslin2 --noZeroInflate --RandomEffect > ${methodName}_noZeroInflate_RandomEffect_${metadataType}_${nSubjects}_${nPerSubject}_${nMicrobes}_${spikeMicrobes}_${nMetadata}_${spikeMetadata}_${effectSize}_${readDepth}.Rout
fi


