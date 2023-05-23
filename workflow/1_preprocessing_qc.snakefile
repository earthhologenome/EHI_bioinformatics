################################################################################
################################################################################
################################################################################
# EHI snakefile for preprocess/QC
# Raphael Eisenhofer 03/2023
#         .----------------.  .----------------.  .----------------.
#        | .--------------. || .--------------. || .--------------. |
#        | |  _________   | || |  ____  ____  | || |     _____    | |
#        | | |_   ___  |  | || | |_   ||   _| | || |    |_   _|   | |
#        | |   | |_  \_|  | || |   | |__| |   | || |      | |     | |
#        | |   |  _|  _   | || |   |  __  |   | || |      | |     | |
#        | |  _| |___/ |  | || |  _| |  | |_  | || |     _| |_    | |
#        | | |_________|  | || | |____||____| | || |    |_____|   | |
#        | |              | || |              | || |              | |
#        | '--------------' || '--------------' || '--------------' |
#         '----------------'  '----------------'  '----------------'
################################################################################
################################################################################
################################################################################

### Setup sample inputs and config

configfile: "preprocess_qc_config.yaml"

import glob
import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_prb_input.py' script, which pulls the information from AirTable and saves-
## it as 'prb_input.tsv'.

# Get list of samples (EHI numbers)
with open("prb_input.tsv", "r") as f:
    SAMPLE = [line.strip() for line in f]

print("Detected these samples")
print(SAMPLE)
### Code to scale time needed by raw read file sizes
### Scaling is based on benchmark data for ~280 jobs 3/4/2023 RE
import os

def estimate_time_download(wildcards):
    fs_sample = f"{config['workdir']}/{wildcards.sample}_filesize.txt"
    with open(fs_sample, 'r') as f:
        input_size = int(f.read().strip())
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # Multiply by 2, and set time based on 30 MB/s download speed.
    estimate_time_download = ((input_size_gb * 2.1 ) + 12 ) / 1.25
    return int(estimate_time_download)

def estimate_time_fastp(wildcards):
    r1_path = f"{config['workdir']}/{wildcards.sample}_raw_1.fq.gz"
    r2_path = f"{config['workdir']}/{wildcards.sample}_raw_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # Add scaling (* 2.5 is for the Gbp to .gz compressed filesize scaling -- e.g. 3 Gbp sample ~ 1.5 GBytes) 
    estimate_time_fastp = ((input_size_gb * 2.5 ) + 6) / 2
    return int(estimate_time_fastp)

def estimate_time_mapping(wildcards):
    r1_path = f"{config['workdir']}/{wildcards.sample}_trimmed_1.fq.gz"
    r2_path = f"{config['workdir']}/{wildcards.sample}_trimmed_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    estimate_time_mapping = ((input_size_gb * 2 ) + 2) * 12
    return int(estimate_time_mapping)

def estimate_time_nonpareil(wildcards):
    r1_path = f"{config['workdir']}/{wildcards.sample}_M_1.fq"
    r2_path = f"{config['workdir']}/{wildcards.sample}_M_2.fq"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # N.b. for estimate_time_nonpareil we estimate from uncompressed fq
    estimate_time_nonpareil = (input_size_gb + 2) * 2
    return int(estimate_time_nonpareil)


def estimate_time_singlem(wildcards):
    r1_path = f"{config['workdir']}/{wildcards.sample}_M_1.fq.gz"
    r2_path = f"{config['workdir']}/{wildcards.sample}_M_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    estimate_time_singlem = ((input_size_gb) + 7) * 4.5
    return int(estimate_time_singlem)

################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("/projects/ehi/data/REP/{prb}.tsv",
                prb=config["prb"]
        )


include: os.path.join(config["codedir"], "rules/create_PRB_folder.smk")
include: os.path.join(config["codedir"], "rules/get_filesize_erda.smk")
include: os.path.join(config["codedir"], "rules/download_raw.smk")
include: os.path.join(config["codedir"], "rules/fastp.smk")
include: os.path.join(config["codedir"], "rules/get_host_genome.smk")
include: os.path.join(config["codedir"], "rules/map_to_host.smk")
include: os.path.join(config["codedir"], "rules/nonpareil.smk")
include: os.path.join(config["codedir"], "rules/singlem.smk")
include: os.path.join(config["codedir"], "rules/coverm_host.smk")
include: os.path.join(config["codedir"], "rules/upload_prb.smk")
include: os.path.join(config["codedir"], "rules/prb_summary.smk")


onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[prb]}" | mailx -s "{config[prb]} ERROR" EMAIL_ADD
          """)
