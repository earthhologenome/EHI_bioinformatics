################################################################################
################################################################################
################################################################################
# EHI snakefile for automatically uploading metadata and data to the ENA
# Raphael Eisenhofer 06/2024
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

### Setup MAG inputs and config
configfile: "ena_upload.yaml"

import os
import glob
import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_ena_input.py' script, which pulls the information from AirTable and saves-
## it as 'ehi_numbers.tsv'.

# Get list of samples (EHI numbers)
with open("ehi_numbers.tsv", "r") as f:
    EHI = [line.strip() for line in f]

print("Detected these samples")
print(EHI)


def estimate_time_download(wildcards, attempt):
    return attempt * 120

rule all:
    input:
        os.path.join(
            config["workdir"], 
            "pipeline_finished"
        )

include: os.path.join(config["codedir"], "rules/log_ehs_start.smk")
include: os.path.join(config["codedir"], "rules/download_raw_ena.smk")
#include: os.path.join(config["codedir"], "rules/download_ena_mags.smk")
include: os.path.join(config["codedir"], "rules/get_sample_checklists.smk")
include: os.path.join(config["codedir"], "rules/register_sample_ena_upload_reads.smk")
#include: os.path.join(config["codedir"], "rules/register_mags_ena.smk")
#include: os.path.join(config["codedir"], "rules/aspera_upload_mags.smk")
include: os.path.join(config["codedir"], "rules/log_ehs_done.smk")


onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[ehs]}" | mailx -s "{config[ehs]} ERROR" EMAIL_ADD
          """)