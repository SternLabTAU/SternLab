#!/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#PBS -r y
#PBS -q adis
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH 
#PBS -N primer_ids_hiv
#PBS -l mem=4000mb
#PBS -J 1-194


python /sternadi/nobackup/volume1/genomes_quantify/hiv-gag/code/barcode_aligner.py /sternadi/nobackup/volume1/genomes_quantify/hiv-gag/barcode.fasta /sternadi/nobackup/volume1/genomes_quantify/hiv-gag/splitted_r/HIV-GAG.R2.$PBS_ARRAY_INDEX /sternadi/nobackup/volume1/genomes_quantify/hiv-gag/quantify_output/HIV-GAG.R2.$PBS_ARRAY_INDEX.primer_ids
