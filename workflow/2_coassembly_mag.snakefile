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
    (row["PR_batch"], row["EHI_number"], row["Assembly_code"]) for _, row in df.iterrows()
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
include: os.path.join(config["codedir"], "rules/metawrap_binning_coassembly.smk")
include: os.path.join(config["codedir"], "rules/metawrap_refinement_coassembly.smk")
include: os.path.join(config["codedir"], "rules/coverm_coassembly.smk")
include: os.path.join(config["codedir"], "rules/gtdbtk_coassembly.smk")
include: os.path.join(config["codedir"], "rules/coassembly_summary.smk")
include: os.path.join(config["codedir"], "rules/log_coASB_finish.smk")

onerror:
    shell("""
            echo "/projects/ehi/data/RUN/{config[abb]}" | mailx -s "{config[abb]} ERROR" EMAIL_ADD
          """)