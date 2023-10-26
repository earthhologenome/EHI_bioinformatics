################################################################################
################################################################################
################################################################################
# EHI snakefile for coassembly/binning
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


import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_eha_input.py' script, which pulls the information from AirTable and saves-
## it as 'asb_input.tsv'.

df = pd.read_csv("asb_input.tsv", sep="\t")

# Use set to create a list of valid combinations of wildcards.
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

def estimate_time_download(wildcards, attempt):
    return attempt * 20

def estimate_time_assembly(wildcards, attempt):
    assembly_size = pd.read_csv("assembly_size.csv")
    assembly_size_value = assembly_size['assembly_size'].values[0]

    if assembly_size_value > 100:
        return attempt * 1440
    else:
        return attempt * 480

def estimate_time_mapping(wildcards, attempt):
    row = get_row(wildcards)
    metagenomic_bases = row["metagenomic_bases"]
    singlem_fraction = row["singlem_fraction"]
    diversity = row["diversity"]
    C = row["C"]
    gbp_post_mapping = calculate_input_size_gb(row["metagenomic_bases"])
    estimate_time_mapping = -59.69 + (2.14 * diversity) - (0.016 * singlem_fraction) + (4.34 * gbp_post_mapping) + (22.22 * C)
    estimate_time_mapping = max(estimate_time_mapping, 20)
    return attempt * int(estimate_time_mapping)

def estimate_time_binning(wildcards, attempt):
    return attempt * 250

def estimate_time_refinement(wildcards, attempt):
    return attempt * 450

def estimate_time_gtdb(wildcards, attempt):
    return attempt * 300

def estimate_time_upload_bam(wildcards, attempt):
    return attempt * 60

def estimate_time_summary(wildcards, attempt):
    return attempt * 60

def estimate_time_coverm(wildcards, attempt):
    return attempt * 10

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
                config["workdir"], 
                "{combo[2]}_QUAST"
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
include: os.path.join(config["codedir"], "rules/coassembly.smk")
include: os.path.join(config["codedir"], "rules/QUAST_coassembly.smk")
include: os.path.join(config["codedir"], "rules/index_coassembly.smk")
include: os.path.join(config["codedir"], "rules/coassembly_mapping.smk")
include: os.path.join(config["codedir"], "rules/upload_coasb_bam.smk")
include: os.path.join(config["codedir"], "rules/concoct_coassembly.smk")
include: os.path.join(config["codedir"], "rules/metabat2_coassembly.smk")
include: os.path.join(config["codedir"], "rules/maxbin2_coassembly.smk")
include: os.path.join(config["codedir"], "rules/metawrap_refinement_coassembly.smk")
include: os.path.join(config["codedir"], "rules/coverm_coassembly.smk")
include: os.path.join(config["codedir"], "rules/gtdbtk_coassembly.smk")
include: os.path.join(config["codedir"], "rules/coassembly_summary.smk")
include: os.path.join(config["codedir"], "rules/log_coASB_finish.smk")

onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[abb]}" | mailx -s "{config[abb]} ERROR" EMAIL_ADD
          """)