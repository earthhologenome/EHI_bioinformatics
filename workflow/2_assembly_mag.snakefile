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


## Set up dynamic times for rules based on input data:
## values are derived from benchmarks (gbp / time required) see URL TO GITHUB MD *********
def get_row(wildcards):
    return df[
        (df["PR_batch"] == wildcards.PRB) &
        (df["EHI_number"] == wildcards.EHI)
    ].iloc[0]

def calculate_input_size_gb(metagenomic_bases):
    # convert from bytes to gigabytes
    return metagenomic_bases / (1024 * 1024 * 1024)

## Rule-specific time estimations
## This also includes retries by attempt
def estimate_time_download(wildcards):
    row = get_row(wildcards)
    input_size_gb = calculate_input_size_gb(row["metagenomic_bases"])
    estimate_time_download = input_size_gb / 1.4
    return attempt * int(estimate_time_download)

def estimate_time_assembly(wildcards):
    row = get_row(wildcards)
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    nonpareil_estimated_coverage = row["nonpareil_estimated_coverage"]
    input_size_gb = calculate_input_size_gb(metagenomic_bases)
    estimate_time_assembly = -217.6446692 + (13.6639860 * diversity) - (23.1672386 * singlem_fraction) + (11.5142151 * input_size_gb) - (0.5745956 * nonpareil_estimated_coverage)
    return attempt * int(estimate_time_assembly)

# def estimate_time_assembly_attempt(wildcards, attempt):
#     return attempt * estimate_time_assembly(wildcards)

# def estimate_time_assembly_mapping(wildcards):
#     row = get_row(wildcards)
#     input_size_gb = calculate_input_size_gb(row["metagenomic_bases"])
#     estimate_time_assembly_mapping = input_size_gb / 0.13
#     return int(estimate_time_assembly_mapping)

# def estimate_time_assembly_mapping_attempt(wildcards, attempt):
#     return attempt * estimate_time_assembly_mapping(wildcards)

# def estimate_time_binning(wildcards):

# def estimate_time_binning_attempt(wildcards, attempt):
#     return attempt * estimate_time_binning(wildcards)

# def estimate_time_coverm(wildcards):

# def estimate_time_coverm_attempt(wildcards, attempt):
#     return attempt * estimate_time_coverm(wildcards)

# def estimate_time_index_assembly(wildcards):

# def estimate_time_index_assembly_attempt(wildcards, attempt):
#     return attempt * estimate_time_index_assembly(wildcards)

# def estimate_time_gtdb(wildcards):

# def estimate_time_gtdb_attempt(wildcards, attempt):
#     return attempt * estimate_time_gtdb(wildcards)

# def estimate_time_refinement(wildcards):

# def estimate_time_refinement_attempt(wildcards, attempt):
#     return attempt * estimate_time_refinement(wildcards)

# def estimate_time_upload_bam(wildcards):

# def estimate_time_upload_bam_attempt(wildcards, attempt):
#     return attempt * estimate_time_upload_bam(wildcards)



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
#         expand(
#             os.path.join(
#                 config["workdir"], "{combo[0]}_{combo[1]}_{combo[2]}_uploaded"
#             ),
#             combo=valid_combinations,
#         ),
#         expand(
#             os.path.join(config["workdir"], "{abb}_pipeline_finished"),
#             abb=config["abb"],
#         )


include: os.path.join(config["codedir"], "rules/create_ASB_folder.smk")
include: os.path.join(config["codedir"], "rules/download_preprocessed.smk")
include: os.path.join(config["codedir"], "rules/individual_assembly.smk")
include: os.path.join(config["codedir"], "rules/QUAST.smk")
# include: os.path.join(config["codedir"], "rules/index_assembly.smk")
# include: os.path.join(config["codedir"], "rules/assembly_mapping.smk")
# include: os.path.join(config["codedir"], "rules/upload_asb_bam.smk")
# include: os.path.join(config["codedir"], "rules/metawrap_binning.smk")
# include: os.path.join(config["codedir"], "rules/metawrap_refinement.smk")
# include: os.path.join(config["codedir"], "rules/coverm_assembly.smk")
# include: os.path.join(config["codedir"], "rules/gtdbtk.smk")
# include: os.path.join(config["codedir"], "rules/assembly_summary.smk")
# include: os.path.join(config["codedir"], "rules/log_ASB_finish.smk")



##For logging errors
onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[abb]}" | mailx -s "{config[abb]} ERROR" EMAIL_ADD
          """)


