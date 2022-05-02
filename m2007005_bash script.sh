#!/bin/bash #

#Step1: Data organization
# First, need to be in the home directory to create a storage space for my planned work
# To create the main working folder for the workflow called ngs_course, use the command mkdir
mkdir ngs_course
# To create a directory for my DNA-seq analysis part
mkdir ngs_course/dnaseq
# To change to the dnaseq directory, use the command cd
cd ngs_course/dnaseq
#To make four different directories to save and organize my different types of data
#The data folder will be used to store any raw data files; the meta folder will be used to contain any info that describes the samples i have used; the results folder will contain all output files of my workflow; th elog files is to track the commands i am running and the specific parameters i am using. 
mkdir data meta results log
# To verify that I have now created all the required direcotries use the command ls
ls -lF
cd ~/ngs_course/dnaseq/data
# My DNA-seq data will be split into trimmed and untrimmed sequences, therefore to create files for this purpose
mkdir untrimmed_fastq
mkdir trimmed_fastq

#Step 2: Downlading raw files
# Now, remain within the data folder, and to download my raw-fastq data
# will download the data using the wget command which allows non-interactive downloads for free utility from the web.
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#now download the annotation file needed for the alignment step
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# To move all the raw untrimmed_fastq files and bed file from the data directory into the untrimmed_fastq directory
mv *fastq.qz ~/ngs_course/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/ngs_course/dnaseq/data

#Now download the hg19 reference genome sequences from ucsc database using the wget command
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#then move the reference genome file into the data folder
mv hg19.fa.gz ~/ngs_course/dnaseq/data/

#Step 3: Downloading all required dependencies
#Installing the tools required to run my NGS pipeline, return to home directory first
cd ~/
#after I am in the home directory, first is to download the miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./ Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib

#Step 4: Redirection
# Now it is the step to capture our FASTQ files through redirection process, where I will redirect the output from my terminal to a print file so i can look at it later.
# For my redirection step, i will utilize the > command to redirect, and the zcat command to help me display my FASTQ files in a compressed form without having to uncompress it
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
ls -lart

#the zcat command will allow to display the content of the fastq files in a compressed form without uncompressing the file
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq
#To assess the quality metrics of my untrimmed fastq data, I will perform a fastqc step on all of my fastq files.
fastqc NGS0001.R1.fastq NGS0001.R2.fastq

mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

#Step 5: Trimmomatic
#To change  my current directory to the untrimmed fastq location

cd ~/ngs_course/dnaseq/data/untrimmed_fastq

#Now we need to trim the adapter sequences and filter out poor quality score reads in our fastq files, the trimmomatic tool will remove also sequences that fall below the minimum required sequence length
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq \
  -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R \
 ILLUMINACLIP:/home/ubuntu//miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

#Step 6: Quality control of trimmed data
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P NGS0001.R1.fastq_trimmed_R_2P
#Now I have trimmed-quality reads which i can take to the alignment step which is aligning the reads to the reference genome.
#For this I need to make a folder called reference where I can move the reference genome sequence reads file into.
mkdir -p ~/ngs_course/dnaseq/data/reference
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/
#Before aligning, I need to make a reference index using the bwa command which will generate many index files
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
#Use the ls command to check that the re required index files have were generated
ls ~/ngs_course/dnaseq/data/reference

#Now since I created trimmed files from my untrimmed_fastq, I can delete the untrimmed fastq files to empty space in the terminal
rm -r ~/ngs_course/dnaseq/data/untrimmed_fastq

#Step 7: Alignment
#first, make a new directory called aligned data to store the data after alignment is complete
mkdir ~/ngs_course/dnaseq/data/aligned_data

#To perform the alignment step, will use the BWA software package to map the dna sequences against the hg19 reference genome.
# The bwa mem algorithm  is designed for illumina sequence reads from 70bp up to megabases, it is the fastest algorithm and the most accurate.
#Specially for reads 70-100bp, it has the best performance compared to other algorithms.
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

#To make sure I am in the aligned-data folder
cd ~/ngs_course/dnaseq/data/aligned_data

#Now I need to convert all sam files into bam, then sort these files and generate an index. all by using the samtools
samtools view -h -b NGS0001.sam > NGS0001.bam


samtools sort NGS0001.bam > NGS0001_sorted.bam


samtools index NGS0001_sorted.bam

#Now again, stay in the aligned-data folder.
cd ~/ngs_course/dnaseq/data/aligned_data

#the flagsstat counts the number of alignments for each FLAG type.It provides counts for each of 13 read categories depending on bit flags in the FLAG field.
samtools flagstat NGS0001_sorted.bam
#To reports alignment summary statistics we can use the samtools idxstats
samtools idxstats NGS0001_sorted.bam

cd ~/ngs_course/dnaseq/data/aligned_data
#Since we have converted the NGS data into bam files, we can remove the sam files to provide more space for the terminal to run  markduplicates and variant calling
rm NGS0001.sam

#To mark the dupliate reads, can use the picard tool as it examines aligned records in the supplied SAM or BAM dataset to locate duplicate molecules.
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam
#After performing markduplicates, Filter BAM based on mapping quality and bitwise flags using samtools. And Set Minimum MAPQ quality score : 20.
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

#Step 8: Variant calling and filtering
##first we need to decompress the hg19 reference genome file using the zcat
zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa
#Then, i have converted the reference from text format and indexed it with the samtools faidx
#the freebayes is a variant calling detector that is designed to find small SNPs, indels..etc, it requires fasta sequence files and BAM annotation file.
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf
#after freebayes runs, the file is huge, therefore we can compress it to make a VCF file
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

#Now, filtering the VCF files generated by freebayes, this VCF files contains a large amount of bad calls and we can use freebayes information fields to filer this such as QUAL, AO, SAF, RPL, RPR.
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

cd ~/ngs_course/dnaseq/data/aligned_data

##to screen the overlaps between my vcf file and the data of annotated bed file, can use the bedtools intersect
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf -b ../annotation.bed > ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

#Now need to compress the vcf file using the bgzip as indexing with tabix requires a compressed input
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

##index the vcf file with tabix ( tabix is a generic indexer)
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz
#I do not need now the trimmed fastq tiles as I have completed the variant calling step and can empty space on the terminal for the variant annotation step.
rm -r ~/ngs_course/dnaseq/data/trimmed_fastq

#Step 9: annotation
#The result of the variant calling is a vcf files that contains the location and data but it does not contain information on the functional significance of each variant. We can use the ANNOVAR free tool for annotation.
cd ~/
#After, I need to upload the annovar.latest.tar file onto openstack then unpack it and set it up by the tar command. 
tar -zxvf annovar.latest.tar.gz
cd annovar
#before i can use annovar , i need to download all the necessary databses for annotation for example the clinvar data base. 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/

#Now convert the VCF to annovar input format
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_R.avinput
#Now run the annovar table function to get a csv output which i can extract as an excel file via filezilla and this will allow me to observe and study each annotated variant. 
./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_R.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

