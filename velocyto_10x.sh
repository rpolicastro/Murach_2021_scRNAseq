#!/bin/bash

##################
## RNA Velocity ##
##################

cd $PBS_O_WORKDIR

module load singuality/3.5.2

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif \
velocyto run10x \
  -@ $NCORES \
  --samtools-memory 5000 \
  -m genome/new/mm10_rmsk.gtf \
  aligned/new/$SAMPLE \
  genome/new/filtered_custom_Mus_musculus.GRCm38.99.chr/genes/genes.gtf
