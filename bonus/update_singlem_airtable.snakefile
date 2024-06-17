################################################################################
################################################################################
################################################################################
# EHI snakefile updating SingleM values in AirTable for new SingleM versions or-
# metapackages
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

configfile: "singlem_update.yaml"

import os
import glob
import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_read_input_singlem.py' script, which pulls the information from AirTable and saves-
## it as 'read_input.tsv'.

df = pd.read_csv('singlem_input.csv')

EHI = df.iloc[:, 0].tolist()

print("Detected these samples")
print(EHI)

def estimate_time_singlem(wildcards, attempt):
    return attempt * 120


rule all:
    input:
        os.path.join(
            config["workdir"],
            "singlem_updated"
        )

include: os.path.join(config["codedir"], "rules/download_reads_singlem.smk")
include: os.path.join(config["codedir"], "rules/singlem_update.smk")
include: os.path.join(config["codedir"], "rules/patch_singlem_updates.smk")

