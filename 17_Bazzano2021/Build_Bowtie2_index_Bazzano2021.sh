#!/bin/bash

GENOME=/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Resources/Bazzano2021/Reference/ILV1-L.fasta
OUTPUT_BASE=Bazzano2021_ILV1_L_idx

bowtie2-build $GENOME $OUTPUT_BASE