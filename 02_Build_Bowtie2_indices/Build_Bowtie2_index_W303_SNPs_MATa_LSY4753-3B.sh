#!/bin/bash

GENOME=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/S_cerevisiae_reference_genomes/S288C_R64-2-1_W303_SNPs_MATa_hml_hmr_LSY4753-3B.fasta
OUTPUT_BASE=S288C_R64-2-1_W303_SNPs_MATa_hml_hmr_LSY4753-3B_idx

bowtie2-build $GENOME $OUTPUT_BASE