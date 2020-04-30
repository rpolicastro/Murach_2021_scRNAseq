#!/bin/bash

ncores=8

##################
## RNA Velocity ##
##################

cd $PBS_O_WORKDIR

## Prepare Singularity Container
## ----------

module load singuality/3.5.2

if [ ! -f scrnaseq_software_velocyto_0.17.17.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:velocyto_0.17.17
fi

## Run RNA-velocity
## ----------

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif \
velocyto run10x \
  -@ $ncores \
  --samtools-memory 5000 \
  -m genome/new/mm10_rmsk.gtf \
  aligned/new/KY_mononuclear \
  genome/new/filtered_custom_Mus_musculus.GRCm38.99.chr/genes/genes.gtf
