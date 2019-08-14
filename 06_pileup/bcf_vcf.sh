#!/bin/bash
#SBATCH --job-name=bcf_vcf
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

module load bcftools/1.6

bcftools call -mv pileup.out > var.raw.vcf
bcftools filter -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > var.flt.vcf



echo "end of pileup"
date




