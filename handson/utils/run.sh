#! /bin/bash
xhost + 127.0.0.1
DIR=/home/fpozoc/Desktop/docker-workshop
docker run --rm -e DISPLAY=host.docker.internal:0 -v "$DIR/DATA:/DATA" -v "$DIR/RESULTS:/RESULTS/report" -it osvaldogc/ufv:2.0 /bin/bash


## docker cp ID:/RESULTS/s_cerevisiae_chr/ # Docker id

# Cuidado con los permisos

# 3. Pipeline
cd ~
## Run FastQC
fastqc /DATA/s_cerevisiae_chrX_read1.fastq.gz -o /RESULTS -t 3
fastqc /DATA/s_cerevisiae_chrX_read2.fastq.gz -o /RESULTS -t 3

## Run Trimmomatic
printf ">TruSeq_Adapter_Index_11\nGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG" > /RESULTS/illumina_adapter.fasta
java -jar /SOFTWARE/trimmomatic/bin/trimmomatic.jar PE -phred33 /DATA/s_cerevisiae_chrX_read1.fastq.gz /DATA/s_cerevisiae_chrX_read2.fastq.gz /RESULTS/s_cerevisiae_chrX_read1_paired.fastq /RESULTS/s_cerevisiae_chrX_read1_unpaired.fastq /RESULTS/s_cerevisiae_chrX_read2_paired.fastq /RESULTS/s_cerevisiae_chrX_read2_unpaired.fastq ILLUMINACLIP:illumina_adapter.fasta:2:33:20:2:true LEADING:36 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:3
fastqc /RESULTS/s_cerevisiae_chrX_read1_paired.fastq -o /RESULTS -t 3
fastqc /RESULTS/s_cerevisiae_chrX_read2_paired.fastq -o /RESULTS -t 3
echo "Trimmomatic has finished." 

# Bowtie2
echo "Bowtie2 has started..."
bowtie2-build /DATA/s_cerevisiae_chrX.fasta /RESULTS/s_cerevisiae_chrX
export BOWTIE2_INDEXES=/RESULTS
bowtie2 -x /RESULTS/s_cerevisiae_chrX -p 2 -1 /RESULTS/s_cerevisiae_chrX_read1_paired.fastq -2 /RESULTS/s_cerevisiae_chrX_read2_paired.fastq -S /RESULTS/s_cerevisiae_chrX.sam
echo "Bowtie2 has finished."

# Samtools
echo "Samtools has started..."
samtools view -Su /RESULTS/s_cerevisiae_chrX.sam -o /RESULTS/s_cerevisiae_chrX.bam
samtools sort /RESULTS/s_cerevisiae_chrX.bam -o /RESULTS/s_cerevisiae_chrX_sorted.bam # Improving performance of analysis
echo "Samtools has finished."

# Stringtie 
echo "Stringtie has started..."
REFERENCE_ANNOTATION=/DATA/s_cerevisiae_chrX.gff
BAM_FILE=/RESULTS/s_cerevisiae_chrX_sorted.bam
OUTPUT_ALL=/RESULTS/stringtie_output/stringtie_all_output.gtf
OUTPUT_ALL=/RESULTS/stringtie_output/stringtie_all_output.gtf
OUTPUT_COV=/RESULTS/stringtie_output/stringtie_cov_output.gtf
stringtie $BAM_FILE -G $REFERENCE_ANNOTATION -o $OUTPUT_ALL -A $OUTPUT_COV
echo "Stringtie has finished."

echo "The analysis has been completed succesfully"
