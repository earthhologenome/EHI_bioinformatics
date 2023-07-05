################################################################################
################################################################################
################################################################################
# EHI snakefile for MAG annotation
# Raphael Eisenhofer 04/2023
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

configfile: "dereplication_mapping.yaml"

import os
import glob
import pandas as pd

## Setup the list of MAGs that we want, this is pulled from airtable using the
## 'get_derep_mags_airtable.py' script.

df = pd.read_csv("dereped_mags.csv", sep=",")

MAG = list(df.iloc[:, 1])

valid_combinations = set(
    (row["ehm"], row["mag_name"]) for _, row in df.iterrows()
)


print("Detected the following MAGs:")
print(MAG)

rule all:
    input:
        expand(
            os.path.join(
                config["magdir"], 
                "{MAG}_anno.tsv.gz"
            ),
            combo=valid_combinations
        ),    
        os.path.join(
            config["magdir"],
            "MAGs_uploaded"
        )

include: os.path.join(config["codedir"], "rules/download_dereplicated.smk")
include: os.path.join(config["codedir"], "rules/dram.smk")
include: os.path.join(config["codedir"], "rules/upload_dram.smk")