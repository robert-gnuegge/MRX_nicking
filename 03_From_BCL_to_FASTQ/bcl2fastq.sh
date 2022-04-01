#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive FASTQ files from BCL files (BCL -> FASTQ)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date: 01/20/22
# version: 2.0

RUNFOLDER="/media/robert/Elements/Deep_sequencing_data/21-02-08-RG_3/21-02-17-All_run_files/"

OUTFOLDER="/media/robert/Elements/Deep_sequencing_data/21-02-08-RG_3/21-02-17-bcl2fastq_output/"

nohup bcl2fastq  --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --auto-set-to-zero-barcode-mismatches --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 -R $RUNFOLDER -o $OUTFOLDER --interop-dir $OUTFOLDER"InterOp" --no-lane-splitting --use-bases-mask Y63N*,I*,Y12N* --sample-sheet "SampleSheet_no_adpt_trim.csv"

# --mask-short-adapter-reads 0 and --minimum-trimmed-read-length 0 suppress masking of short reads as NNN reads
# --auto-set-to-zero-barcode-mismatches avoids index collision, as some indices differ only at a single base
# --use-bases-mask Y63N*,I*,Y12N* selects how many cycles to analyze for first, index, and second read
# --sample-sheet $OUTFOLDER"SampleSheet_no_adpt_trim.csv" specifies the modified sample sheet, where adapter definitions were removed, to suppress adapter trimming