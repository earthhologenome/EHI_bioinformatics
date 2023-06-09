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
    (row["PR_batch"], row["EHI_number"], row["Assembly_code"], row["metagenomic_bases"], row["singlem_fraction"], row["diversity"], row["C"]) for _, row in df.iterrows()
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
def estimate_time_download(wildcards, attempt):
    row = get_row(wildcards)
    input_size_gb = calculate_input_size_gb(row["metagenomic_bases"])
    estimate_time_download = (input_size_gb / 1.4)
    estimate_time_download = max(estimate_time_download, 2)
    return attempt * int(estimate_time_download)

def estimate_time_assembly(wildcards, attempt):
    row = get_row(wildcards)
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    C = row["C"]
    gbp_post_mapping = calculate_input_size_gb(metagenomic_bases)
    estimate_time_assembly = -112.01 + (11.10 * diversity) - (51.77 * singlem_fraction) + (12.72 * gbp_post_mapping) - (53.05 * C)
    estimate_time_assembly = max(estimate_time_assembly, 20)
    return attempt * int(estimate_time_assembly)

def estimate_time_quast(wildcards, attempt):
    return attempt * 5

def estimate_time_assembly_mapping(wildcards, attempt):
    row = get_row(wildcards)
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    C = row["C"]
    gbp_post_mapping = calculate_input_size_gb(row["metagenomic_bases"])
    estimate_time_assembly_mapping = -59.69 + (2.14 * diversity) - (0.016 * singlem_fraction) + (4.34 * gbp_post_mapping) + (22.22 * C)
    estimate_time_assembly_mapping = max(estimate_time_assembly_mapping, 10)
    return attempt * int(estimate_time_assembly_mapping)

def estimate_time_binning(wildcards, attempt):
    row = get_row(wildcards)
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    nonpareil_estimated_coverage = row["nonpareil_estimated_coverage"]
    C = row["C"]
    gbp_post_mapping = calculate_input_size_gb(row["metagenomic_bases"])
    estimate_time_binning = -595.79 + (25.59 * diversity) - (0.18 * singlem_fraction) + (5.72 * gbp_post_mapping) + (169.88 * C)
    estimate_time_binning = max(estimate_time_binning, 20)
    return attempt * int(estimate_time_binning)

def estimate_time_refinement(wildcards, attempt):
    return attempt * 210

def estimate_time_gtdb(wildcards, attempt):
    return attempt * 90

def estimate_time_upload_bam(wildcards, attempt):
    return attempt * 5

def estimate_time_summary(wildcards, attempt):
    return attempt * 60


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



##For logging errors
onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[abb]}" | mailx -s "{config[abb]} ERROR" EMAIL_ADD
          """)


