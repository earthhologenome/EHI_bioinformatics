################################################################################
################################################################################
################################################################################
# EHI snakefile for MAG dereplication and count table creation
# Raphael Eisenhofer 05/2023
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

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_read_input.py' script, which pulls the information from AirTable and saves-
## it as 'read_input.tsv'.

df = pd.read_csv("read_input.tsv", sep="\t")

# Use set to create a list of valid combinations of wildcards. Note that 'ID' = EHA number.
valid_combinations = set(
    (row["PR_batch"], row["EHI_number"]) for _, row in df.iterrows()
)

## Setup the list of MAGs that we want, this is pulled from airtable using the
## 'get_mags_airtable.py' script.

df2 = pd.read_csv("mags.csv", sep=",")

MAG = list(df2.iloc[:, 0])

# Remove the first element of the list (which is the column name "genome")
MAG.pop(0)

print("Detected the following MAGs:")
print(MAG)


rule all:
    input:
        os.path.join(
            config["workdir"],
            "tables_uploaded"
        )

include: os.path.join(config["codedir"], "rules/create_DMB_folder.smk")
include: os.path.join(config["codedir"], "rules/download_preprocessed.smk")
include: os.path.join(config["codedir"], "rules/download_mags.smk")
include: os.path.join(config["codedir"], "rules/drep.smk")
include: os.path.join(config["codedir"], "rules/gtdbtk_full_tree.smk")
include: os.path.join(config["codedir"], "rules/index_mags.smk")
include: os.path.join(config["codedir"], "rules/mag_mapping.smk")
include: os.path.join(config["codedir"], "rules/coverm_mags.smk")
include: os.path.join(config["codedir"], "rules/upload_tables.smk")
