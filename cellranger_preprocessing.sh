#!/bin/bash

############################
## Cell Ranger Processing ##
############################

cd $PBS_O_WORKDIR

## Prepare Singularity
## ----------

module load singularity/3.5.2

if [ ! -f scrnaseq_software_cellranger_3.1.0.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:cellranger_3.1.0
fi

## Run Cell Ranger
## ----------

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_cellranger_3.1.0.sif \
cellranger count \
    --id KY_mononuclear \
    --fastqs sequences \
    --sample KY_mononuclear \
    --transcriptome genome/new/filtered_custom_Mus_musculus.GRCm38.99.chr \
    --expect-cells 10000

mv KY_mononuclear aligned/new
