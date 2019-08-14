#!/bin/bash
#SBATCH --job-name=pileup
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load samtools/1.7

samtools mpileup -uf ../02_reference/chr22.fa ../05_alignment/trimmed_NA12878_sort.bam > pileup.out


