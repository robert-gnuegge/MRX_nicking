#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive deduplicated alignments for NGS reads (FASTQ -> BAM)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/27/22
# version: 2.0

SAMPLE=$1
READDIR="/media/robert/Elements/Deep_sequencing_data/21-07-28-RG_4/Merged_FastQ_files"
IDX="S288C_R64-2-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr_PHO5-SrfIcs_idx"
export BOWTIE2_INDEXES=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bowtie2_indices
ADAPTERS=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/2nd_Adapter_linkers/2nd_Adpt_linkers.fa

# To keep time
SECONDS=0
Elapsed () {
	TIME="[$(printf "%02d" $(($SECONDS /3600))):$(printf "%02d" $(($SECONDS / 60))):$(printf "%02d" $(($SECONDS % 60)))]"
}

Elapsed
printf "\n$TIME Read2 length filtering..."
fastp --in1 "$READDIR"/"$SAMPLE"/"$SAMPLE"_R1.fastq.gz --in2 "$READDIR"/"$SAMPLE"/"$SAMPLE"_R2.fastq.gz \
--out1 tmp_R1.fastq.gz --out2 tmp_R2.fastq.gz --length_required=12 \
--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g \
-h fastp_1.html 2>> fastp_1.log

Elapsed
printf "\n$TIME UMI extraction..."
fastp --in1 tmp_R1.fastq.gz --in2 tmp_R2.fastq.gz --out1 tmp2_R1.fastq.gz --out2 tmp2_R2.fastq.gz \
--umi --umi_loc=read2 --umi_len=12 \
--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering \
-h fastp_2.html 2>> fastp_2.log

rm tmp_R1.fastq.gz tmp_R2.fastq.gz tmp2_R2.fastq.gz

Elapsed
printf "\n$TIME Preprocessing and mapping Read1; converting and sorting alignments..."
fastp --in1 tmp2_R1.fastq.gz --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_fasta $ADAPTERS -h fastp_3.html 2>>fastp_3.log --stdout \
| bowtie2 --local --very-sensitive -p 4 -x $IDX - 2> Bowtie2.log \
| samtools view -bu - \
| samtools sort - -o tmp.bam

rm tmp2_R1.fastq.gz

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

rm tmp.bam
rm tmp.bam.bai

Elapsed
printf "\n$TIME Indexing..."
samtools index "$SAMPLE".bam

Elapsed
printf "\n$TIME Done.\n"

# Notes -------------------------------------------------------------------
# * Read2 has to be length-filtered, as truncated UMIs break umi_tools dedup.
# * UMIs are moved from Read2 to the Read1 headers. Subsequently, only Read1 is processed further.
# * To evaluate the effect of using UMIs for deduplication, we also use umi_tools dedup with the --ignore-umi parameter. This can be skipped.