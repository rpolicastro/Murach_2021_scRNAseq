#!/bin/bash

cd $PBS_O_WORKDIR
cd ..

module load singularity

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_cellranger_3.1.0.sif \
  cellranger aggr \
  --id=4day \
  --csv=scripts/aggr_csv/4day.csv

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_cellranger_3.1.0.sif \
  cellranger aggr \
  --id=14day \
  --csv=scripts/aggr_csv/14day.csv
