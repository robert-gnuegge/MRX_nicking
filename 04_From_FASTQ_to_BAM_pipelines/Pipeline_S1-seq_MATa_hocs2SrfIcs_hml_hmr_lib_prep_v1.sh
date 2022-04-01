#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive deduplicated alignments for NGS reads (FASTQ -> BAM)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 12/20/21
# version: 2.0

SAMPLE=$1
READDIR="/media/robert/Elements/Deep_sequencing_data/19-05-21-RG_2/19-05-21-FastQ_files"
READS="$READDIR"/"$SAMPLE".fastq.gz
IDX="S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr_idx"
export BOWTIE2_INDEXES=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bowtie2_indices

# To keep time
SECONDS=0
Elapsed () {
	TIME="[$(printf "%02d" $(($SECONDS /3600))):$(printf "%02d" $(($SECONDS / 60))):$(printf "%02d" $(($SECONDS % 60)))]"
}

Elapsed
printf "\n$TIME Preprocessing, converting, sorting, and mapping..."
fastp --in1 $READS --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--umi --umi_loc=read1 --umi_len=12 --umi_skip=10 2>>fastp.log --stdout \
| bowtie2 --local --very-sensitive -p 4 -x $IDX - 2> Bowtie2.log \
| samtools view -bu - \
| samtools sort - -o tmp.bam

Elapsed
printf "\n$TIME Indexing..."
samtools index tmp.bam

# <can be skipped>
# to evaluate UMI effect
Elapsed
printf "\n$TIME Deduplicating without UMIs..."
umi_tools dedup  --stdin=tmp.bam --umi-separator=":" --no-sort-output --ignore-umi --log=umi_tools_dedup_ignore-umi.log --stdout=tmp2.bam
rm tmp2.bam
# </can be skipped>

Elapsed
printf "\n$TIME Deduplicating and sorting..."
umi_tools dedup  --stdin=tmp.bam --umi-separator=":" --no-sort-output --log=umi_tools.log \
| samtools sort - -o "$SAMPLE".bam

Elapsed
printf "\n$TIME Indexing..."
samtools index "$SAMPLE".bam

# Elapsed
printf "\n$TIME Removing temporary files..."
rm tmp.bam
rm tmp.bam.bai

Elapsed
printf "\n$TIME Done.\n"

# Notes -------------------------------------------------------------------
# * Processing the input files in one pipe is not possible here, as umi_tools dedup cannot read from stdin. That is why we create a temporary file (tmp.bam) which is later removed.
# * To evaluate the effect of using UMIs for deduplication, we also use umi_tools dedup with the --ignore-umi parameter. This can be skipped.