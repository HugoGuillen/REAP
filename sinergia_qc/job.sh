#!/bin/bash
#$ -l h_cpu=12:00:00
#$ -l h_vmem=1G
#$ -cwd
#$ -o $JOB_ID_output.txt
#$ -e $JOB_ID_error.txt
#$ -N snakemake
#$ -M carlos.pulido@dbmr.unibe.ch
#$ -m ea
#$ -q all.q
#$ -pe smp 1

snakemake --printshellcmds --cluster "qsub -q all.q -pe smp 8 -o /data/projects/p283_rna_and_disease/projects/sinergia/sinergia/qc/output -e /data/projects/p283_rna_and_disease/projects/sinergia/sinergia/qc/error -N QC" --jobs 10 --latency-wait 3600 
