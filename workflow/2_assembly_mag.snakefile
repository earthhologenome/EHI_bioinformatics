################################################################################
################################################################################
################################################################################
# EHI snakefile for assembly/binning and MAG annotation (individual assemblies)
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
    (row["PR_batch"], row["EHI_number"], row["Assembly_code"], row["metagenomic_bases"], row["singlem_fraction"], row["diversity"], row["nonpareil_estimated_coverage"]) for _, row in df.iterrows()
)


################################################################################
### Setup the desired outputs
rule all:
    input:
        os.path.join(
                config["workdir"], 
                "ERDA_folder_created"
        ),
        expand(
            os.path.join(
                config["workdir"], "{combo[0]}_{combo[1]}_{combo[2]}_QUAST"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"], "{combo[0]}_{combo[1]}_{combo[2]}_uploaded"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(config["workdir"], "{abb}_pipeline_finished"),
            abb=config["abb"],
        )


include: os.path.join(config["codedir"], "rules/create_ASB_folder.smk")
include: os.path.join(config["codedir"], "rules/download_preprocessed.smk")
include: os.path.join(config["codedir"], "rules/individual_assembly.smk")
include: os.path.join(config["codedir"], "rules/QUAST.smk")
include: os.path.join(config["codedir"], "rules/index_assembly.smk")
include: os.path.join(config["codedir"], "rules/assembly_mapping.smk")
include: os.path.join(config["codedir"], "rules/upload_asb_bam.smk")
include: os.path.join(config["codedir"], "rules/metawrap_binning.smk")
include: os.path.join(config["codedir"], "rules/metawrap_refinement.smk")
include: os.path.join(config["codedir"], "rules/coverm_assembly.smk")
include: os.path.join(config["codedir"], "rules/gtdbtk.smk")
include: os.path.join(config["codedir"], "rules/assembly_summary.smk")
include: os.path.join(config["codedir"], "rules/log_ASB_finish.smk")


## Set up dynamic times for rules based on input data:

### Define the dynamic time estimates based on input file sizes
## values are derived from benchmarks (gbp / time required) see URL TO GITHUB MD *********
def estimate_time_download(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_download = input_size_gb / 1.4
    return int(estimate_time_download)

def estimate_time_assembly(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    nonpareil_estimated_coverage = row["nonpareil_estimated_coverage"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_assembly = -217.6446692 + (13.6639860 * diversity) - (23.1672386 * singlem_fraction) + (11.5142151 * input_size_gb) - (0.5745956 * nonpareil_estimated_coverage)
    return int(estimate_time_assembly)

def estimate_time_assembly_attempt(wildcards, attempt)
    return attempt * estimate_time_assembly

def estimate_time_assembly_mapping(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_assembly_mapping = input_size_gb / 0.13
    return int(estimate_time_assembly_mapping)

def estimate_time_binning(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_bining = input_size_gb / 0.03
    return int(estimate_time_bining)

def estimate_time_coverm(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_coverm = input_size_gb / 2.5
    return int(estimate_time_coverm)

def estimate_time_index_assembly(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_index_assembly = input_size_gb / 0.9
    return int(estimate_time_index_assembly)

def estimate_time_gtdb(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_gtdb = input_size_gb / 0.03
    return int(estimate_time_gtdb)

def estimate_time_refinement(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_refinement = input_size_gb / 0.015
    return int(estimate_time_refinement)

def estimate_time_upload_bam(wildcards):
    row = df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]
    metagenomic_bases = row["metagenomic_bases"]
    # convert from bytes to gigabytes
    input_size_gb = metagenomic_bases / (1024 * 1024 * 1024)
    estimate_time_upload_bam = input_size_gb / 2
    return int(estimate_time_upload_bam)



##For logging errors
onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[abb]}" | mailx -s "{config[abb]} ERROR" EMAIL_ADD
          """)


