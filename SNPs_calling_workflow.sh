#!  usr/bin/bash
#Begin Variant/SNPs calling workflow"
#Make sure you are in your home directory to start workflow  ..........."

mkdir ~/Documents/SRR
echo "Where are the sequence Reads stored? ....." 
read sequence
cp  $sequence/SRR* ~/Documents/SRR

echo "Creating Working Directory"

echo "Unzipping Raw fastq files from Sequence Read"
for i in *.zip
   do
   unzip ${i}
done

cd ~
echo "Creating Working Directory"
cd  ~/Documents/SRR

echo "Copying Raw Reads in Downloads"
cp ~/Downloads/SRR* .

echo "Unzipping Raw fastq files from Sequence Read"
for i in *.zip
   do
   unzip ${i}
done

mkdir -p ~/Documents/SRR/data/untrimmed_fastq
find . -name *.fastq.gz -exec cp -t ./data {} +
cd ~/Documents/SRR/data
mv *.fastq.gz  ~/Documents/SRR/data/untrimmed_fastq
 cd ~/Documents/SRR/data/untrimmed_fastq

echo "Running FastQC ..."
fastqc *.fastq*

mkdir -p ~/Documents/SRR/results/fastqc_untrimmed_reads

echo "Saving FastQC results..."
mv *.zip ~/Documents/SRR/results/fastqc_untrimmed_reads/
mv *.html ~/Documents/SRR/results/fastqc_untrimmed_reads/

cd ~/Documents/SRR/results/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done

echo "Saving summary..."
cat */summary.txt > ~/Documents/SRR/fastqc_summaries.txt

cd ~/Documents/SRR/data/untrimmed_fastq

echo "Running Trimmomatic ..."
cp ~/Desktop/NexteraPE-PE.fa  ~/Documents/SRR/data/untrimmed_fastq

for infile in *_1.fastq.gz
   do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done

mkdir  ~/Documents/SRR/trimmed_fastq
cd   ~/Documents/SRR/data/untrimmed_fastq
 mv *.trim* ~/Documents/SRR/trimmed_fastq
 cd ~/
 
# Download your reference genome into your result directory
 curl -L -o  Desktop/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz 
 gunzip ~/Desktop/ref_genome/ecoli_rel606.fasta.gz
cd ~/Documents/SRR/results/

genome=~/Desktop/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

cd ~/Documents/SRR/trimmed_fastq
mv *un.trim* ~/Documents/SRR/ 

echo "Unzipping Raw fastq files in trimmed_fastq"
for i in *.gz
   do
   gunzip ${i}
done

for fq1 in ~/Documents/SRR/trimmed_fastq/*_1.trim.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.fastq)
    echo "base name is $base"

    fq1=~/Documents/SRR/trimmed_fastq/${base}_1.trim.fastq
    fq2=~/Documents/SRR/trimmed_fastq/${base}_2.trim.fastq
    sam=~/Documents/SRR/results/sam/${base}.aligned.sam
    bam=~/Documents/SRR/results/bam/${base}.aligned.bam
    sorted_bam=~/Documents/SRR/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/Documents/SRR/results/bcf/${base}_raw.bcf
    variants=~/Documents/SRR/results/vcf/${base}_variants.vcf
    final_variants=~/Documents/SRR/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
done

for i in ~/Documents/SRR/results/vcf/*_final_variants.vcf
    do
    echo ${i}
    grep -v "#" ${i} | wc -l
done

echo "Your work is done"

