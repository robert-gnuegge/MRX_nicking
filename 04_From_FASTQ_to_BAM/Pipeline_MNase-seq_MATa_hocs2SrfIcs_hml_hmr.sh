#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive deduplicated alignments for NGS reads (FASTQ -> BAM)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 11/07/21
# version: 1.0

# File locations
READDIR=/media/robert/Elements/Deep_sequencing_data/21-11-04-RG_7/Merged_FastQ_files

# Sample name
SAMPLE=$1
printf "\n========== Processing $SAMPLE =========="

READ1="$READDIR"/"$SAMPLE"/"$SAMPLE"_R1.fastq.gz
READ2="$READDIR"/"$SAMPLE"/"$SAMPLE"_R2.fastq.gz
export BOWTIE2_INDEXES=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/Bowtie2_indices

# To keep run time
SECONDS=0
Elapsed () {
	ELAPSEDTIME="$(printf "%02d" $(($SECONDS /3600))):$(printf "%02d" $(($SECONDS / 60))):$(printf "%02d" $(($SECONDS % 60)))"
}

# Process reads, map, and deduplicate
printf "\n[$(date '+%T')] Starting analysis..."
printf "\n[$(date '+%T')] Piping fastp -> bowtie2 -> samtools view -> samtools collate"
printf "\n-> sammtools fixmate -> samtools sort -> samtools markdup..."

fastp --in1 $READ1 --in2 $READ2 --stdout --detect_adapter_for_pe \
	--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	2> fastp.log \
	| bowtie2 --very-sensitive --no-mixed --no-discordant -p 4 \
	-x $IDX \
	--interleaved - \
	2> Bowtie2.log \
	| samtools view -bu - \
	| samtools collate -Ou - \
	| samtools fixmate -m - - \
	| samtools sort - \
	| samtools markdup -rs 2> Samtools_markdup.log - "$SAMPLE".bam

# Index BAM file
printf "\n[$(date '+%T')] Indexing final BAM file..."
samtools index "$SAMPLE".bam

# end
Elapsed
printf  "\n[$(date '+%T')] Done (run time: $ELAPSEDTIME).\n"

# Notes -------------------------------------------------------------------
# * samtools view -b results in BAM output format
# * samtools collate -O outputs to stdout
# * samtools fixmate -m adds mate score tags, which samtools markdup needs to pick best alignments to keep
# * samtools view/collate -u prevents BAM output suppression to save time
# * samtools markdup -r removes duplicates, instead of marking them