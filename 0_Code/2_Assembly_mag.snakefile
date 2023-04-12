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
            os.path.join(
                config["workdir"], "{combo[0]}", "{combo[1]}", "{combo[2]}_QUAST"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}",
                "{combo[1]}",
                "{combo[2]}_assembly_coverM.txt",
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}",
                "{combo[1]}",
                "{combo[2]}",
                "gtdbtk/classify/gtdbtk.bac120.summary.tsv"
            ),
            combo=valid_combinations,
        )


include: os.path.join(config["codedir"], "rules/create_ASB_folder.smk")
include: os.path.join(config["codedir"], "rules/download_preprocessed.smk")
include: os.path.join(config["codedir"], "rules/individual_assembly.smk")
include: os.path.join(config["codedir"], "rules/QUAST.smk")
include: os.path.join(config["codedir"], "rules/index_assembly.smk")
include: os.path.join(config["codedir"], "rules/assembly_mapping.smk")
include: os.path.join(config["codedir"], "rules/metawrap_binning.smk")
include: os.path.join(config["codedir"], "rules/metawrap_refinement.smk")
include: os.path.join(config["codedir"], "rules/coverm_assembly.smk")
include: os.path.join(config["codedir"], "rules/gtdbtk.smk")
include: os.path.join(config["codedir"], "rules/dram.smk")


################################################################################
### Generate output summary table
# rule generate_summary:
#     input:
#         coverm = "{config['workdir']}/{PRB}/{EHI}/{EHA}_assembly_coverM.txt",
#     output:
#         "3_Outputs/{group}_coassembly_summary.tsv"
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     params:
#         group = "{group}"
#     threads:
#         1
#     resources:
#         mem_gb=16,
#         time='00:05:00'
#     message:
#         "Creating final coassembly summary table"
#     shell:
#         """
#         #Create the final output summary table
#         #parse QUAST outputs for assembly stats
#         echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\tassembly_mapping_percent" > headers.tsv


#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do cat 3_Outputs/2_Coassemblies/{params.group}_QUAST/{params.group}_assembly_report.tsv >> {params.group}_temp_report.tsv;
#         done


#         #Add in the % mapping to assembly stats
#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do echo $(basename ${{sample/.bam/}}) >> {params.group}_sample_ids.tsv;
#         done

#         paste {params.group}_sample_ids.tsv {params.group}_temp_report.tsv > {params.group}_temp2_report.tsv

#         #Add in the # of bins
#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do cat {params.group}_bins.tsv >> {params.group}_number_bins.tsv;
#         done

#         paste {params.group}_temp2_report.tsv {params.group}_number_bins.tsv > {params.group}_temp3_report.tsv

#         ls -l 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam | wc -l > {params.group}_n_samples.tsv

#         nsamples=$( cat {params.group}_n_samples.tsv )
#         nsamples1=$(( nsamples + 1 ))
#         echo $nsamples1
#         for sample in `seq 2 $nsamples1`;
#             do cut -f"$sample" 3_Outputs/6_Coassembly_CoverM/{params.group}_assembly_coverM.txt | sed -n 3p >> {params.group}_relabun.tsv;
#         done

#         paste {params.group}_temp3_report.tsv {params.group}_relabun.tsv > {params.group}_temp4_report.tsv
#         #Combine them into the final assembly report
#         cat headers.tsv {params.group}_temp4_report.tsv > {output}
#         """
# ################################################################################
# ### Clean up
# rule clean:
#     input:
#         expand("3_Outputs/{group}_coassembly_summary.tsv", group=GROUPS.keys())
#     output:
#         "3_Outputs/pipeline_complete.txt"
#     threads:
#         1
#     resources:
#         mem_gb=16,
#         time='00:05:00'
#     message:
#         "Cleaning up temp files"
#     shell:
#         """
#         rm *.tsv
#         touch {output}
#         """
