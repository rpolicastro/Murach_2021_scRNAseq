#!/bin/bash

## Variables

samples=(KY_mononuclear, pax7_dta_4day, pax7_tdt_4day, tdt_parental)

##########################################
## Submit Cell Ranger Processing Script ##
##########################################

## Prepare Singularity
## ----------

module load singularity/3.5.2

if [ ! -f scrnaseq_software_cellranger_3.1.0.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:cellranger_3.1.0
fi

## Submit Jobs
## ----------

for sample in ${samples[@]}; do
  qsub -v SAMPLE=${sample} -l nodes=1:ppn=8,vmem=128gb,walltime=24:00:00 cellranger_preprocessing.sh
done
