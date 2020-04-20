#!/bin/bash

##################
## RNA Velocity ##
##################

## Prepare Singularity Container
## ----------

module load singuality

if [ ! -f scrnaseq_software_velocyto_exec.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:velocyto_exec
fi

## Run RNA-velocity
## ----------


