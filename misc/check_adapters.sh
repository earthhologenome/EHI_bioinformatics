#!/bin/sh
#SBATCH -c 2 --mem 32G # number of cores

#This is a simple slurm script for testing whether you're using the right adapters in fastp.
#It calls fastp, then fastqc on the output. You can then view the fastqc .html and check for adapter contamination there.

#Load conda environment
#N.B., you'll need to create a conda environment first that contains fastp and fastqc
source activate fastp

#Run fastp
        fastp \
            --in1 *_1.fastq.gz --in2 *_2.fastq.gz \
            --out1 test_trimmed_1.fastq.gz --out2 test_trimmed_2.fastq.gz \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread 2 \
            --html test.fastp_html \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

#Run fastqc
fastqc test_trimmed*.fastq.gz -p 2
