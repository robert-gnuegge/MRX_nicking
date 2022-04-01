#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive deduplicated alignments for NGS reads (FASTQ -> BAM)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 02/13/22
# version: 1.0

# File locations
READS=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Bazzano2021/Reads/SRR13769867_exo1_sgs1_1.fastq.gz
LEADER=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Bazzano2021/Leader/ILV1-L_leader_sequence.fa
BARCODES=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Bazzano2021/Barcodes/SRR13769867_Barcodes.fasta
IDX=Bazzano2021_ILV1_L_idx
export BOWTIE2_INDEXES=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bowtie2_indices

# To keep time
SECONDS=0
Elapsed () {
	TIME="[$(printf "%02d" $(($SECONDS /3600))):$(printf "%02d" $(($SECONDS / 60))):$(printf "%02d" $(($SECONDS % 60)))]"
}

Elapsed
printf "\n$TIME Starting analysis..."

Elapsed
printf "\n$TIME Removing leader..."
flexbar -r $READS --barcodes $LEADER --barcode-error-rate 0.1 --threads 4  --target tmp --zip-output GZ

mv tmp.log Flexbar_leader_removal.log

Elapsed
printf "\n$TIME Extracting UMIs and barcodes..."
umi_tools extract --stdin=tmp_barcode_leader.fastq.gz --bc-pattern=NNNNNNNNNNNN --log=umi_tools_extract.log \
| flexbar -r - --barcodes $BARCODES --zip-output GZ --target tmp --threads 4

mv tmp.log Flexbar_debarcoding.log
rm tmp_barcode_leader.fastq.gz
find . -name 'tmp_barcode*' | xargs -I{} rename 's/tmp_barcode/exo1_sgs1/' {}

declare -a arr=("0" "35" "60" "75" "90") 
for TIME in "${arr[@]}"
do
  
  SAMPLE=exo1_sgs1_T"$TIME"
  
  Elapsed
  printf "\n\n$TIME Processing $SAMPLE..."
  
  Elapsed
  printf "\n$TIME Mapping, converting, and sorting..."
  bowtie2 -U $SAMPLE.fastq.gz --very-sensitive -p 4 -x $IDX 2> Bowtie2_$SAMPLE.log \
  | samtools view -bu - \
  | samtools sort - -o tmp.bam
  
  Elapsed
  printf "\n$TIME Indexing..."
  samtools index tmp.bam
  
  Elapsed
  printf "\n$TIME Deduplicating and sorting..."
  umi_tools dedup --stdin=tmp.bam --no-sort-output --log=umi_tools_dedup_$SAMPLE.log --method=unique \
  | samtools sort - -o "$SAMPLE".bam

  rm tmp.bam
  rm tmp.bam.bai
  
  Elapsed
  printf "\n$TIME Indexing..."
  samtools index "$SAMPLE".bam

done

Elapsed
printf "\n$TIME Done.\n"

# Notes:
# * use of umi_tools dedup with option --method=unique is suboptimal, as this does not consider PCR and sequencing errors,
#   however, as there are so many reads at the same position (e.g. uncut), using a different method takes too long to run
#   (see https://github.com/CGATOxford/UMI-tools/issues/241 for a similar case)