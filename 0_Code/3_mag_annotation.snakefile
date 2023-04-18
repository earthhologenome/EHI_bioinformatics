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

### Setup sample inputs and config

configfile: "assembly_mag_config.yaml"

import glob
import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_eha_input.py' script, which pulls the information from AirTable and saves-
## it as 'asb_input.tsv'.

df = pd.read_csv("asb_input.tsv", sep="\t")

# Use set to create a list of valid combinations of wildcards. Note that 'ID' = EHA number.
valid_combinations = set(
    (row["PR_batch"], row["EHI_number"], row["ID"]) for _, row in df.iterrows()
)


rule all:
    input:
        # expand(
        #     os.path.join(
        #         config["workdir"], 
        #         "{combo[0]}/", 
        #         "{combo[1]}/", 
        #         "{combo[2]}/", 
        #         "DRAM/", 
        #         "{MAG}_anno.tsv.gz"
        #     ),
        #     combo=valid_combinations
        # ),
        expand(os.path.join(
                config["workdir"],
                "{combo[0]}/", 
                "{combo[1]}/", 
                "{combo[2]}/", 
                "DRAM/",
                "MAGs_uploaded"
            ),
            combo=valid_combinations
        )

include: os.path.join(config["codedir"], "rules/dram.smk")
include: os.path.join(config["codedir"], "rules/upload_mags.smk")