#!/bin/bash

############################
## Cell Ranger Processing ##
############################

cd $PBS_O_WORKDIR

## Run Cell Ranger
## ----------

module load singularity/3.5.2

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_cellranger_3.1.0.sif \
cellranger count \
    --id $SAMPLE \
    --fastqs sequences/$SAMPLE \
    --sample $SAMPLE \
    --transcriptome genome/new/filtered_custom_Mus_musculus.GRCm38.99.chr \
    --expect-cells 10000

#mv $sample aligned/new
