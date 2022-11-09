# EHI_bioinformatics
Bioinformatics pipeline to process EHI data

### General information:
This pipeline uses snakemake, and manages dependencies using conda (or mamba) for reproducibility and deployability. The 0_Code directory contains the snakefiles, scripts, and conda environment yamls. 

## Getting started:
Firstly, you'll want to clone this directory to the system where you want to run the analyses:
```
git clone https://github.com/anttonalberdi/EHI_bioinformatics.git
```

## Preprocessing pipeline (1_Preprocess_QC.snakefile)
*updated 14/10/2022, Raphael Eisenhofer*

### What this pipeline does:
This step of the pipeline quality filters (including adapter trimming, polyX tail removal) metagenomic reads using fastp (0.23.1). Host genomes are indexed using BowTie2 (2.4.4), before the fastp-filtered reads being mapped to the host genome/s using BowTie2 (default settings). Nonpareil (3.4.1) is then run on the host-filtered reads to estimate metagenome diversity and assembly coverage. CoverM (0.6.1) is run on the host-mapped BAMs to calculate the read counts against the host genome/s. Finally, the summary statistics of the preprocessing and host mapping are collated into a final report.

Here is a simplified DAG (directed acyclic graph) of the above steps:

![1_Preprocess_QC](figures/dag_1_Preprocess_QC.png)

### Usage:
Currently, the snakefile searches for .fastq.gz files located in this path (assuming you are launching the snakefile from the current directory):
```
2_Reads/1_Untrimmed/*_1.fastq.gz
```
There are a couple of options for getting your data here:
- 1) Create symbolic links. This means you don't have to move or copy the files:
`ln -s reads_folder/*.fastq.gz 2_Reads/1_Untrimmed/`
- 2) You can just put the reads in this directory.

(note that the fastq file suffixes should be **'_1.fastq.gz'** and **'_2.fastq.gz'**).

Next, you'll need to do the same thing for your host reference genome/s, symbolic linking them or placing them in this folder:
```
1_References/
```
Note that the host reference genome/s need to be gzip compressed -- **'.gz'** suffix (default when downloading from NCBI). 

That's all the setup required to get the pipeline running. Now you just need to launch the snakefile using snakemake. How you do this depends on your HPC server job queueing system. For Mjolnir, I use the following:
```
snakemake \
-s 0_Code/1_Preprocess_QC.snakefile \
-j 10 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} --cores {threads}" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600
```

`snakemake -s 0_Code/1_Preprocess_QC.snakefile -j 10 --cluster "sbatch --mem {resources.mem_gb}G --cores {threads} --time {resources.time} {threads}" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 600`

I recommend adding the `--dry-run` or `-n` command to the above code initially, as this will let you figure out if everything is working as expected.

I've written the pipeline such that it handles the requesting of optimised resources (RAM/CPUs) for each job based on the specific snakemake rule. The `--conda-prefix` snakemake option tells snakemake to look for conda environment in a particular directory, which saves having to reinstall them each time you run from a new directory. 

Here's a illustrative summary of each rule and it's input files and output files:

![1_Preprocess_QC](figures/file_structure_1_Preprocess_QC.png)
