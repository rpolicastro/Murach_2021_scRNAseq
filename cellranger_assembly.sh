#!/bin/bash

ncores=8

############################
## Custom Genome Assembly ##
############################

if [ ! -d genome/new ]; then
  mkdir -p genome/new
fi

## Prepare Files
## ----------

cd genome/new

## GTF annotation.

if [ ! -f Mus_musculus.GRCm38.99.chr.gtf ]; then
  wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.chr.gtf.gz
  gunzip Mus_musculus.GRCm38.99.chr.gtf.gz
fi

if [ ! -f custom_Mus_musculus.GRCm38.99.chr.gtf ]; then
  cat Mus_musculus.GRCm38.99.chr.gtf Ai9_tdTomato.gtf > custom_Mus_musculus.GRCm38.99.chr.gtf
fi

## FASTA assembly.

if [ ! -f Mus_musculus.GRCm38.dna_sm.primary_assembly.fa ]; then
  wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
  gunzip Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
fi

if [ ! -f custom_Mus_musculus.GRCm38.dna_sm.primary_assembly.fa ]; then
  cat Mus_musculus.GRCm38.dna_sm.primary_assembly.fa Ai9_tdTomato.fasta > custom_Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
fi

## Prepare Singularity Container
## ----------

cd ../..

module load singularity/3.5.2

if [ ! -f scrnaseq_software_cellranger_3.1.0.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:cellranger_3.1.0
fi

## Assemble Genome
## ----------

## Keep only protein coding genes.

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_cellranger_3.1.0.sif \
cellranger mkgtf \
  genome/new/custom_Mus_musculus.GRCm38.99.chr.gtf \
  genome/new/filtered_custom_Mus_musculus.GRCm38.99.chr.gtf \
  --attribute=gene_biotype:protein_coding

## Create the reference genome.

cd genome/new

singularity exec -eCB `pwd` -H `pwd` ../../scrnaseq_software_cellranger_3.1.0.sif \
cellranger mkref \
  --genome=filtered_custom_Mus_musculus.GRCm38.99.chr \
  --fasta=custom_Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
  --genes=filtered_custom_Mus_musculus.GRCm38.99.chr.gtf \
  --nthreads=$ncores \
  --memgb=64

cd ../..
