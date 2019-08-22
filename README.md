# Variant Calling with Samtools (Basics)  
   
  
This repository is a usable, publicly available tutorial for introduction to basics of variant calling. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.    
   
   
  
## Overview
This tutorial is designed to introduce the tools, data types and workflow of variant detection. We will align reads to the genome, look for differences between reads and reference genome sequence, and filter the detected genomic variation manually to understand the computational basis of variant calling.

We cover the concepts of detecting small variants (SNVs and indels) in human genomic DNA using a small set of reads from chromosome 22.

### Learning Objectives
At the end of the course, you will be able to:

Work with the FASTQ format and base quality scores
Align reads to generate a BAM file and subsequently generate a pileup file
Run the FreeBayes variant caller to find SNVs and indels
Visualise BAM files using the Integrative Genomics Viewer (IGV) and identify likely SNVs and indels by eye


### Note on Source Data
The tutorial is based on analysis of short read data from the exome of chromosome 22 of a single human individual. This individual, a Caucasian female from the US, has been sequenced many times and is one of the best genetically characterised humans in the world - her family was one of the original 30 trios in Hapmap project (www.hapmap.org), she was part of the 1000 Genomes project (www.1000genomes.org), the 10Gen project (www.sequenceontology.org/resources/10Gen.html),  and is often used for benchmarking NGS tools (eg http://www.completegenomics.com/sequence-data/download-data/, http://www.broadinstitute.org/gsa/wiki/index.php/NA12878_test_data).    

Contents   
1.    [Data transfer](#data-transfer)
2.    [Inspecting FASTQ files](#inspecting-fastq-files)
3.    [Quality evaluation using FASTQC](#quality-evaluation-using-fastqc)
4.    [Quality control using sickle](#quality-control-using-sickle)
5.    [Aligning reads to a genome using bwa](#aligning-reads-to-a-genome-using-bwa)
6.    [Visualise the BAM file with IGV](#visualise-the-bam-file-with-igv)
7.    [Generate a pileup file](#generate-a-pileup-file)
8.    [Filtering using BCFtools](#filtering-using-bcftools)





<h2 id="First_Point_Header">Data transfer</h2>

Lets make few directories for the tutorial
```bash
./Variant-Calling-with-Samtools/
    ├── 01_data/
    ├── 02_reference/
```

Then change the directory to the **Variant-Calling-with-Samtools/** using:  
```bash
cd Variant-Calling-with-Samtools 
```  

Copy the data and reference from the source location, if it has not been already done for you.

```bash
cp /UCHC/PublicShare/VariantWorkshop/data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq ./01_data/
cp /UCHC/PublicShare/VariantWorkshop/reference/chr22*  ./02_reference/

```     

Once the files have been copied to the above two folders the folder structure will look like:  
```
01_data/
└── NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
02_reference/
├── chr22.amb
├── chr22.ann
├── chr22.bwt
├── chr22.fa
├── chr22.pac
└── chr22.sa
```  

 


<h2 id="First_Point_Header">Inspecting FASTQ files</h2>


This set of commands will allow you to open the FASTQ file and inspect the data.  Each read is represented by 4 lines.  Details on the FASTQ file format can be found <a href="https://bio Informatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>.   
 

Change the directory to the `data/` folder using:   
```bash
cd 01_data/  
```  

To inspect the first few lines in the FASTQ file can use the `head` command, as follows which will print the first few lines into the terminal window.     
<pre style="color: silver; background: black;">-bash-4.2$ head NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
@61CC3AAXX100125:7:118:2538:5577/1
GACACCTTTAATGTCTGAAAAGAGACATTCACCATCTATTCTCTTGGAGGGCTACCACCTAAGAGCCTTCATCCCC
+
?>CADFEEEBEDIEHHIDGGGEEEEHFFGIGIIFFIIEFHIIIIHIIFFIIIDEIIGIIIEHFFFIIEHIFA@?==
@61CC3AAXX100125:7:1:17320:13701/1
CTCAGAAGACCCTGAGAACATGTGCCCAAGGTGGTCACAGTGCATCTTAGTTTTGTACATTTTAGGGAGATATGAG
+
?BCAAADBBGGHGIDDDGHFEIFIIIIFGEIFIIFIGIGEFIIGGIIHEFFHHHIHEIFGHHIEFIIEECE?>@89
@61CC3AAXX100125:7:93:5100:14497/1
CTCAACTGGCTGAAAGTATTATCAATAGAAAGGAATGTTCAGGTTCTTCAATTTTAGAGTGCCCTGGCCTAGAAGA
</pre>  

 
The command below will let you inspect more part of the file.
```
-bash-4.2$ less NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
```   

 
The command below will help you count number of reads in the file
``` 
-bash-4.2$ grep -c "^@61CC" NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
```   


<h2 id="Second_Point_Header">Quality evaluation using FASTQC</h2>

Lets create a directory called **fastqc** folder inside the main folder **var-intro**, if it has not created for you. Then go into the fastqc folder using changing directory.   
```bash
mkdir 03_fastqc
```   
So now the folder structure will look like:  
```
Variant-Calling-with-Samtools/
├── 01_data
├── 02_reference
└── 03_fastqc
```   

change the folder to the **03_fastqc/** 
```bash
cd 03_fastqc/
```   

In order to evalaute the general quality of reads in the file we will be using FASTQC package.  The command can take multiple files as input and outputs a quality report in html format. Once the file is generated you have to transfer the file to your local computer to open it and examine the results carefully.   
  
```bash
#!/bin/bash 
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load fastqc 
fastqc -o . ../01_data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
```  
The above full script can also be found in here which is called [fastqc.sh](/03_fastqc/fastqc.sh).   

This will produce a fastqc report in HTML format which need to downloaded to your local computer to view. To copy the file from the Xanadu cluster please use the `transfer.cam.uchc.edu` node.   

```bash
scp user_name@transfer.cam.uchc.edu:/FULL_PATH_to_FILE/NA12878.GAIIx.exome_chr22.1E6reads.76bp_fastqc.html . 
```  

The report will have a summary of your FASTQ file.
![](/images/fastqc_summary.png)   

Then it will consists of a report on the quality of reads as;   
![](/images/fastqc_quality_score.png)   
 
Also it will consists of other information which you are welcome to explor during the workshop.  


 
  
<h2 id="Third_Point_Header">Quality control using sickle</h2>

Lets make a directory called *04_trimmed_reads* in the main folder, and then change the working directory to that folder. Now the folder structure will look like:  
```
Variant-Calling-with-Samtools/
├── 01_data
├── 02_reference
├── 03_fastqc
└── 04_trimmed_reads
```   

Sickle performs quality control on illumina paired-end and single-end short read data using a sliding window. As the window slides along the fastq file, the average score of all the reads contained in the window is calculated. Should the average window score fall beneath a set threshold, <a href="https://github.com/najoshi/sickle/blob/master/README.md">sickle</a> determines the reads responsible and removes them from the run. After visiting the SRA pages for our data, we see that our data are single end reads. Let's find out what sickle can do with these:

<pre style="color: silver; background: black;">-bash-4.2$ module load sickle

-bash-4.2$ sickle

<strong>Usage</strong>: sickle <command> [options]

<strong>Command</strong>:
pe      paired-end sequence trimming
se      single-end sequence trimming

--help, display this help and exit
--version, output version  Information and exit</pre>

We have single-end sequences. 

<pre style="color: silver; background: black;">-bash-4.2$ sickle se

Usage: sickle se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>

Options:
-f, --fastq-file, Input fastq file (required)
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)
-o, --output-file, Output trimmed fastq file (required)
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-x, --no-fiveprime, Don't do five prime trimming.
-n, --trunc-n, Truncate sequences at position of first N.
-g, --gzip-output, Output gzipped files.
--quiet, Don't print out any trimming information
--help, display this help and exit
--version, output version information and exit</pre>   
  
The quality may be any score from 0 to 40. The default of 20 is much too low for a robust analysis. We want to select only reads with a quality of 35 or better. Additionally, the desired length of each read is 75bp. Again, we see that a default of 20 is much too low for analysis confidence. We want to select only reads whose lengths exceed 45bp. 

Once we have performed data trimming we will recheck the quality of data using fastqc.
 

Let's put all of this together for our sickle script using our downloaded fastq files:

```
-bash-4.2$ nano sickle_run.sh
```
                                                                                             
```
#!/bin/bash 
#SBATCH --job-name=sickle_run
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load sickle

sickle se -t sanger -f ../01_data/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq \ 
	-o trimmed_NA12878.fastq \
	-l 45 \
	-q 25 

module unload sickle

module load fastqc
fastqc -t 4 -o  ../03_fastqc/ trimmed_NA12878.fastq

```                                                                                              

The full scrip can be is called [sickle_run.sh](/04_trimmed_reads/sickle_run.sh). To run the script using `sbatch` command type:
```
sbatch sickle_run.sh 
```   

This will result in trimmed reads in the **04_trimmed_reads/** folder:   
```
04_trimmed_reads/
├── sickle_run.sh
└── trimmed_NA12878.fastq
```   

Then by running *fastqc* program it will analyze the trimmed read files and produce the output files in **03_fastqc** folder:   
```
03_fastqc/
├── trimmed_NA12878_fastqc.html
└── trimmed_NA12878_fastqc.zip
```   

 

  
<h2 id="Fourth_Point_Header">Aligning reads to a genome using bwa</h2>   

Lets make a directory called **05_alignment** and will then change the working directory to that folder. So now the folder structure will look like:   
  
```
Variant-Calling-with-Samtools/
├── 01_data
├── 02_reference
├── 03_fastqc
├── 04_trimmed_reads
└── 05_alignment
```    
  

The basic process here to map individual reads - from the input sample FASTQ file - to a matching region on the reference genome. We will execute in two steps

1. Align the reads with BWA
We will map (align) the reads with the BWA tool to the human reference genome. For this tutorial, use Human reference genome 19 (hg19) - this is hg19 from UCSC. The indexes are located are /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bwa .

2. Sort the bam file based on coordinates of alignment


### Assess the alignment data

We can generate some mapping statistics from the BAM file to assess the quality of our alignment. We will use two commands that are part of `samtools`.

Run IdxStats.
IdxStats generates a tab-delimited output with four columns. Each line consists of a reference sequence name (e.g. a chromosome), reference sequence length, number of mapped reads and number of placed but unmapped reads.
We can see that most of the reads are aligning to chromosome 22 as expected.

Run Flagstat. 
Note that in this case the statistics are not very informative. This is because the dataset has been generated for this workshop and much of the noise has been removed (and in fact we just removed a lot more noise in the previous step); also we are using single ended read data rather than paired-end so some of the metrics are not relevant.



```bash
nano bwa_run.sh
```
  
```bash
#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#export TMPDIR=/home/CAM/$USER/tmp/

## STEP1
module load bwa
bwa mem -t 8 ../02_reference/chr22 ../04_trimmed_reads/trimmed_NA12878.fastq -o trimmed_NA12878.sam


## STEP2
module load samtools
samtools view -@ 8 -bhS trimmed_NA12878.sam -o trimmed_NA12878.bam
samtools sort -@ 8 trimmed_NA12878.bam -o trimmed_NA12878_sort.bam
samtools index -b trimmed_NA12878_sort.bam

samtools idxstats trimmed_NA12878_sort.bam >> align_stat.txt
samtools flagstat trimmed_NA12878_sort.bam >> align_stat.txt 

```
  
Command
```
-t : number of processors been used
```

The full slurm scrip is called [bwa_align.sh](/05_alignment/bwa_align.sh) which can be found in the **05_alignment** folder. 
Now you can run this using ` sbatch bwa_run.sh`

Once the command finishes we will inspect the alignment stats by opening `align_stat.txt` file using `less align_stat.txt` command. 
The final file structure will look like: 
```
05_alignment/
├── align_stat.txt
├── bwa_align.sh
├── trimmed_NA12878.bam
├── trimmed_NA12878.sam
├── trimmed_NA12878_sort.bam
└── trimmed_NA12878_sort.bam.bai
```   




<h2 id="Fifth_Point_Header">Visualise the BAM file with IGV</h2>

To visualise the alignment data:

Transfer the bam file and the index file to your local machine, using `transfer.cam.uchc.edu` transfer node. 
```bash
scp user_id@transfer.cam.uchc.edu:/UCHC/PublicShare/Variant_Detection_Tutorials/Variant-Calling-with-Samtools/05_alignment/trimmed_NA12878_sort.bam* .
```  

Start IGV on your computer. You will need Java installed on your computer to run IGV)open the BAM file in your current IGV session.
Once IGV opens, it will show you the BAM file. This may take a bit of time as the data is downloaded.
Our reads for this tutorial are from chromosome 22, so select chr22 from the second drop box under the toolbar. 
Zoom in to view alignments of reads to the reference genome.

![](/images/igv1.jpg)
Try looking at region `chr22:36,006,744-36,007,406`

Can you see a few variants?

![](/images/igv_mb.jpg) 
  
<h2 id="Sixth_Point_Header">Generate a pileup file</h2>

Lets make a directory called **06_pileup/** if its not already made out for you, and then change the directory to that folder, where the folder strcutrue will look like:   
```   
Variant-Calling-with-Samtools/
├── 01_data
├── 02_reference
├── 03_fastqc
├── 04_trimmed_reads
├── 05_alignment
└── 06_pileup
```    


A pileup is essentially a column-wise representation of the aligned reads - at the base level - to the reference. The pileup file summarises all data from the reads at each genomic region that is covered by at least one read. Each row of the pileup file gives similar information to a single vertical column of reads in the IGV view.

The current generation of variant calling tools do not output pileup files, and you don't need to do this section in order to use FreeBayes in the next section. However, a pileup file is a good illustration of the evidence the variant caller is looking at internally, and we will produce one to see this evidence.
We will use `mpileup` function of samtools for creating a pileup output.  The general command usage syntax is

<strong>
`samtools mpileup [-EB] [-C capQcoef] [-r reg] [-f in.fa] [-l list] [-Q minBaseQ] [-q minMapQ] in.bam [in2.bam [...]] ` </strong>

Tip: The pileup file we generated has 10 columns:
 1. chromosome
 2. position
 3. current reference base
 4. consensus base from the mapped reads
 5. consensus quality
 6. SNV quality
 7. maximum mapping quality
 8. coverage
 9. bases within reads
 10. quality values   
  
Further information on (9):
Each character represents one of the following (the longer this string, higher the coverage):

. = match on forward strand for that base
, = match on reverse strand
ACGTN = mismatch on forward
acgtn = mismatch on reverse
+[0-9]+[ACGTNacgtn]+' = insertion between this reference position and the next
-[0-9]+[ACGTNacgtn]+' = deletion between this reference position and the next
^ = start of read
$ = end of read
BaseQualities = one character per base in ReadBases, ASCII encoded Phred score

Following samtools we will use `bcftools` to fiter the variants.


```bash
#!/bin/bash
#SBATCH --job-name=pileup
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date
#export TMPDIR=/home/CAM/$USER/tmp/

module load samtools
samtools mpileup -uf ../02_reference/chr22.fa ../05_alignment/trimmed_NA12878_sort.bam > pileup.out

```
The full slurm scrip is called [pileup.sh](/06_pileup/pileup.sh).   


<h2 id="Seventh_Point_Header">Filtering using BCFtools</h2>
BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. A more detailed instruction can be found here (https://samtools.github.io/bcftools/bcftools.html)

```
-m, --multiallelic-caller  alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c)
-v, --variants-only  output variant sites only
-e exclude
%QUAL<20 || DP>100   : Call SNPs and short INDELs, then mark low quality sites and sites with the read depth exceeding a limit. (The read depth should be adjusted to about twice the average read depth as higher read depths usually indicate problematic regions which are often enriched for artefacts.)
```   

following is the slurm script.  
```
#!/bin/bash
#SBATCH --job-name=bcf_vcf
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools/1.6

bcftools call -mv pileup.out > var.raw.vcf
bcftools filter -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > var.flt.vcf
```

The script is called [bcf_vcf.sh](/06_pileup/bcf_vcf.sh). Once the script has been sucessfully executed using `sbatch` command, following files will be resulted as the output:   
```
06_pileup/
├── bcf_vcf.sh
├── pileup.out
├── pileup.sh
├── var.raw.vcf
└── var.flt.vcf
```   


The `var.flt.vcf` file is the VCF file after the filtering and the first few lines are shown bellow and you are welcome to explore this file using the `less` command as shown bellow:    


```bash
less var.flt.vcf 
```   
 

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ../05_alignment/trimmed_NA12878_sort.bam
chr22   16069707        .       C       G       10.7923 LowQual DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60   GT:PL   1/1:40,3,0
chr22   16075294        .       G       A       46.8798 PASS    DP=9;VDB=0.392388;SGB=-0.511536;RPB=0.809011;MQB=0.924584;MQSB=0.974597;BQB=0.924584;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=3,3,2,1;MQ=60   GT:PL   0/1:80,0,165
chr22   16078283        .       G       A       7.30814 LowQual DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60   GT:PL   1/1:36,3,0
chr22   16079289        .       A       G       5.04598 LowQual DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60   GT:PL   1/1:33,3,0
chr22   16081478        .       T       G       10.7923 LowQual DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60   GT:PL   1/1:40,3,0
chr22   16084621        .       T       C       40.1806 PASS    DP=25;VDB=0.701568;SGB=-0.590765;RPB=0.677057;MQB=1;MQSB=1;BQB=0.16378;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=11,9,1,4;MQ=60        GT:PL   0/1:75,0,255
chr22   16084718        .       G       A       222     PASS    DP=48;VDB=0.699348;SGB=-0.693021;RPB=0.999946;MQB=0.854333;MQSB=0.964251;BQB=0.103423;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=18,3,19,8;MQ=59        GT:PL   0/1:255,0,211
chr22   16084859        .       G       A       88      PASS    DP=33;VDB=0.0758337;SGB=-0.680642;RPB=0.317069;MQB=1;BQB=8.08456e-05;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,21,0,12;MQ=60 GT:PL   0/1:121,0,183
chr22   16087780 
```   



  

