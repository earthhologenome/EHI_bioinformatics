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


import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'get_eha_input.py' script, which pulls the information from AirTable and saves-
## it as 'asb_input.tsv'.

df = pd.read_csv("asb_input.tsv", sep="\t")

# Use set to create a list of valid combinations of wildcards. Note that 'ID' = EHA number.
valid_combinations = set(
    (row["PR_batch"], row["EHI_number"], row["ID"]) for _, row in df.iterrows()
)


################################################################################
### Setup the desired outputs
rule all:
    input:
        expand(
            os.path.join(config["workdir"], "{abb}_ERDA_folder_created"),
            abb=config["abb"],
        ),
        expand(
            os.path.join(config["workdir"], "{combo[0]}", "{combo[1]}_M_1.fq.gz"),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(config["workdir"], "{combo[0]}", "{combo[1]}_M_2.fq.gz"),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"], "{combo[0]}/" "{combo[1]}/" "{combo[2]}_contigs.fasta"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[1]}_{combo[2]}.bam",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"], "{combo[0]}/", "{combo[1]}/", "{combo[1]}_{combo[2]}_uploaded"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"], "{combo[0]}/", "{combo[1]}/", "{combo[2]}_QUAST"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}_binning/binning_complete",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}_refinement/",
                "{combo[2]}_metawrap_50_10_bins.stats",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}_assembly_coverM.txt",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}/",
                "gtdbtk/classify/gtdbtk.bac120.summary.tsv",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"], 
                "{combo[0]}/", 
                "{combo[1]}/", 
                "{combo[2]}/", 
                "DRAM/", 
                "{MAG}_anno.tsv.gz"
            ),
            combo=valid_combinations,
            MAG=range(1, 3000)
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}_final_stats.tsv",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(config["workdir"], "{abb}_pipeline_finished"),
            abb=config["abb"],
        ),


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
include: os.path.join(config["codedir"], "rules/dram.smk")
include: os.path.join(config["codedir"], "rules/assembly_summary.smk")
include: os.path.join(config["codedir"], "rules/log_ASB_finish.smk")
