#! /bin/bash

echo "Creating directories in $PWD..."

DIR=$PWD/pipeline
SOFTWARE=$DIR/software
DATA=$DIR/data
RESULTS=$DIR/results 

mkdir -p $DIR
echo "$DIR created!"

mkdir -p $SOFTWARE 
echo "$SOFTWARE created!"

mkdir -p $DATA
echo "$DATA created!"

mkdir -p $RESULTS
echo "$RESULTS created!"

# 1. Data
echo "Downloading data..."
wget -q -P $DATA/ ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz 
echo "Fasta file downloaded here $DATA"
wget -q -P $DATA/ ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz 
echo "Annotation file downloaded here: $DATA"
gunzip $DATA/*.gz
echo "Files extracted!"
wget -q -P $DATA/ https://github.com/willblev/assembly_workshop_MA_2016/raw/master/s_cerevisiae_chrX_read1.fastq.gz
echo "FastQC-1 ChrX downloaded here: $DATA"
wget -q -P $DATA/ https://github.com/willblev/assembly_workshop_MA_2016/raw/master/s_cerevisiae_chrX_read2.fastq.gz
echo "FastQC-2 ChrX downloaded here: $DATA"
wget -q -P $DATA/ https://github.com/willblev/assembly_workshop_MA_2016/raw/master/s_cerevisiae_chrX.fasta
echo "Fasta file ChrX Downloaded here: $DATA"
wget -q -P $DATA/ https://github.com/willblev/assembly_workshop_MA_2016/raw/master/s_cerevisiae_chrX.gff 
echo "Annotation file ChrX Downloaded here: $DATA"
mv $DATA/GCF_000146045.2_R64_genomic.fna $DATA/s_cerevisiae.fasta
mv $DATA/GCF_000146045.2_R64_genomic.gff $DATA/s_cerevisiae.gff
echo "Data stored here: $DATA"

# 2. Software
## Install FastQC
echo "Downloading FastQC..."
wget -q -P $SOFTWARE/ http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip -q $SOFTWARE/fastqc_v0.11.5.zip -d $SOFTWARE && rm $SOFTWARE/fastqc*.zip
echo "FastQC downloaded here: $SOFTWARE"
FASTQC=$SOFTWARE/FastQC
export PATH=$FASTQC:$PATH
chmod 755 $FASTQC
echo "fastqc executable has been located here: $FASTQC"

## Install Trimmomatic
echo "Downloading Trimmomatic..."
wget -q -P $SOFTWARE/ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
echo "Trimmomatic downloaded here: $SOFTWARE"
unzip -q $SOFTWARE/Trimmomatic-0.36.zip -d $SOFTWARE && rm $SOFTWARE/Trimmomatic-0.36.zip
echo "Trimmomatic is ready to be executed with: java -jar $SOFTWARE/Trimmomatic-0.36/trimmomatic-0.36.jar --help"

## Install Bowtie2
echo "Downloading Bowtie2..."
wget -q -P $SOFTWARE/ https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip
echo "Bowtie2 downloaded here: $SOFTWARE"
unzip -q $SOFTWARE/bowtie2-2.3.3.1-linux-x86_64.zip -d $SOFTWARE && rm $SOFTWARE/bowtie2-2.3.3.1-linux-x86_64.zip
BOWTIE2=$SOFTWARE/bowtie2-2.3.3.1-linux-x86_64
export PATH=$BOWTIE2:$PATH
echo "bowtie2 executable has been located here: $FASTQC"

## Install Samtools
# sudo apt-get install samtools
wget -q -P $SOFTWARE/ https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xjf $SOFTWARE/samtools-1.9.tar.bz2 -C $SOFTWARE && rm $SOFTWARE/samtools-1.9.tar.bz2
cd $SOFTWARE/samtools-1.9/ && ./configure --prefix=$SOFTWARE/samtools-1.9
cd $SOFTWARE/samtools-1.9/ && make
cd $SOFTWARE/samtools-1.9/ && make install
SAMTOOLS=$SOFTWARE/samtools-1.9/bin
export PATH=$SAMTOOLS:$PATH
cd $DIR

## Install Stringtie
echo "Downloading Stringtie..."
wget -q -P $SOFTWARE/ http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.Linux_x86_64.tar.gz
echo "Stringtie downloaded here: $SOFTWARE"
tar xvfz $SOFTWARE/stringtie-1.3.0.Linux_x86_64.tar.gz -C $SOFTWARE && rm $SOFTWARE/stringtie-1.3.0.Linux_x86_64.tar.gz
STRINGTIE=$SOFTWARE/stringtie-1.3.0.Linux_x86_64
export PATH=$STRINGTIE:$PATH
echo "stringtie executable has been located here: $STRINGTIE"

# 3. Pipeline
## Run FastQC
echo "FastQC has started..."
FASTQC_DIR=$RESULTS/fastqc
mkdir -p $FASTQC_DIR
fastqc $DATA/s_cerevisiae_chrX_read1.fastq.gz -o $FASTQC_DIR -t 3
fastqc $DATA/s_cerevisiae_chrX_read2.fastq.gz -o $FASTQC_DIR -t 3
echo "FastQC has finished." 

## Run Trimmomatic
echo "Trimmomatic has started..."
TRIMMOMATIC_DIR=$RESULTS/trimmomatic
mkdir -p $TRIMMOMATIC_DIR
printf ">TruSeq_Adapter_Index_11\nGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG" > $SOFTWARE/Trimmomatic-0.36/illumina_adapter.fasta
java -jar $SOFTWARE/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE -phred33 $DATA/s_cerevisiae_chrX_read1.fastq.gz \
$DATA/s_cerevisiae_chrX_read2.fastq.gz \
$TRIMMOMATIC_DIR/s_cerevisiae_chrX_read1_paired.fastq \
$TRIMMOMATIC_DIR/s_cerevisiae_chrX_read1_unpaired.fastq \
$TRIMMOMATIC_DIR/s_cerevisiae_chrX_read2_paired.fastq \
$TRIMMOMATIC_DIR/s_cerevisiae_chrX_read2_unpaired.fastq \
ILLUMINACLIP:illumina_adapter.fasta:2:33:20:2:true LEADING:36 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:3
fastqc $TRIMMOMATIC_DIR/s_cerevisiae_chrX_read1_paired.fastq -o $FASTQC_DIR -t 3
fastqc $TRIMMOMATIC_DIR/s_cerevisiae_chrX_read2_paired.fastq -o $FASTQC_DIR -t 3
echo "Trimmomatic has finished." 

# Bowtie2
echo "Bowtie2 has started..."
BOWTIE2_DIR=$RESULTS/trimmomatic
mkdir -p $BOWTIE2_DIR
bowtie2-build $DATA/s_cerevisiae_chrX.fasta $BOWTIE2_DIR/s_cerevisiae_chrX
export BOWTIE2_INDEXES=$BOWTIE2_DIR
bowtie2 -x $BOWTIE2_DIR/s_cerevisiae_chrX -p 2 -1 $TRIMMOMATIC_DIR/s_cerevisiae_chrX_read1_paired.fastq -2 $TRIMMOMATIC_DIR/s_cerevisiae_chrX_read2_paired.fastq -S $BOWTIE2_DIR/s_cerevisiae_chrX.sam
echo "Bowtie2 has finished."

# Samtools
echo "Samtools has started..."
SAMTOOLS_DIR=$RESULTS/samtools
mkdir -p $SAMTOOLS
samtools view -Su $BOWTIE2_DIR/s_cerevisiae_chrX.sam -o $SAMTOOLS/s_cerevisiae_chrX.bam
samtools sort $SAMTOOLS/s_cerevisiae_chrX.bam -o $SAMTOOLS/s_cerevisiae_chrX_sorted.bam # Improving performance of analysis
echo "Samtools has finished."

# Stringtie 
echo "Stringtie has started..."
REFERENCE_ANNOTATION=$DATA/s_cerevisiae_chrX.gff
BAM_FILE=$RESULTS/s_cerevisiae_chrX_sorted.bam
OUTPUT_ALL=$RESULTS/stringtie/stringtie_all_output.gtf
OUTPUT_ALL=$RESULTS/stringtie/stringtie_all_output.gtf
OUTPUT_COV=$RESULTS/stringtie/stringtie_cov_output.gtf
stringtie $BAM_FILE -G $REFERENCE_ANNOTATION -o $OUTPUT_ALL -A $OUTPUT_COV
echo "Stringtie has finished."
echo "The analysis has been completed succesfully"
