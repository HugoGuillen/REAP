#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set max wallclock time
#SBATCH --time=24:00:00

# set name of job
#SBATCH --job-name=SM

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=carlos.pulido@dbmr.unibe.ch

# set memory [K|M|G|T]
#SBATCH --mem=1G

# set working directory "--chdir=./"
#SBATCH

# set error file
#SBATCH --error="%A.error"

# set output file
#SBATCH --output="%A.output"


# commands to execute
# source .snakemake_venv/bin/activate
snakemake -s harvest --printshellcmds --cluster "sbatch --mem 90G -N 1 -o /data/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/logs/output -e /data/projects/p283_rna_and_disease/projects/sinergia/sinergia/sinergia_main/logs/error -J SM_" --jobs 10 --latency-wait 3600
