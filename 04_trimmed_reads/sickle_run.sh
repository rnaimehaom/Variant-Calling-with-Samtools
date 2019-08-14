#!/bin/bash
#SBATCH --job-name=sickle_run
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname

module load sickle

sickle pe -t sanger -f ../01_data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq \
        -o trimmed_NA12878.fastq \
        -l 45 \
        -q 25  

module unload sickle

module load fastqc
fastqc -t 4 -o  ../03_fastqc/ trimmed_NA12878.fastq


