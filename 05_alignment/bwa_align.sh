#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

## STEP 1
module load bwa/0.7.17
bwa mem -t 8 ../02_reference/chr22 ../04_trimmed_reads/trimmed_NA12878.fastq -o trimmed_NA12878.sam
module unload bwa/0.7.17


## STEP 2
module load samtools/1.7
samtools view -@ 8 -bhS trimmed_NA12878.sam -o trimmed_NA12878.bam
samtools sort -@ 8 trimmed_NA12878.bam -o trimmed_NA12878_sort.bam
samtools index -b trimmed_NA12878_sort.bam

samtools idxstats trimmed_NA12878_sort.bam >> align_stat.txt
samtools flagstat trimmed_NA12878_sort.bam >> align_stat.txt


