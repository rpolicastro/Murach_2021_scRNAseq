#!/bin/bash

samples=(KY_mononuclear pax7_dta_4day pax7_tdt_4day tdt_parental)
ncores=8

################################
## RNA-velocity Preprocessing ##
################################

## Prepare Sigularity Container
## ----------

module load singularity

if [ ! -f scrnaseq_software_velocyto_0.17.17.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:velocyto_0.17.17
fi

## 10X Chromium Samples
## ----------

for sample in ${samples[@]}; do
  qsub -v SAMPLE=${sample},NCORES=${ncores} -l nodes=1:ppn=8,vmem=128gb,walltime=24:00:00 velocyto_10x.sh
done
