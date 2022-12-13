################################################################################
################################################################################
################################################################################
# EHI snakefile for exploring non-bacterial/archaeal DNA in metagenomes
# Raphael Eisenhofer 12/2022
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

configfile: "0_Code/configs/eukaryotic_assessment_config.yaml"

### Setup sample inputs
import os
from glob import glob

GROUP = [ dir for dir in os.listdir('2_Reads/4_Host_removed')
         if os.path.isdir(os.path.join('2_Reads/4_Host_removed', dir)) ]

SAMPLE = [os.path.relpath(fn, "2_Reads/4_Host_removed").replace("_M_1.fastq.gz", "")
            for group in GROUP
            for fn in glob(f"2_Reads/4_Host_removed/{group}/*_1.fastq.gz")]

print("Detected these sample groups:")
print(GROUP)
print("Detected the following samples:")
print(SAMPLE)

################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("3_Outputs/{group}_eukaryotic_report.tsv", group=GROUP)

################################################################################
### Classify contigs using CAT
rule cat:
    input:
        contigs = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        diamond = temp("3_Outputs/2_Coassemblies/{group}/{group}.alignment.diamond"),
        faa = "3_Outputs/2_Coassemblies/{group}/{group}.predicted_proteins.faa",
        gff = "3_Outputs/2_Coassemblies/{group}/{group}.predicted_proteins.gff",
        classif = "3_Outputs/2_Coassemblies/{group}/{group}.contig2classification.txt",
        final_output = "3_Outputs/2_Coassemblies/{group}/{group}.CAT_final_output.tsv"
    params:
        group = "{group}",
        database = expand("{database}", database=config['database']),
        taxonomy = expand("{taxonomy}", taxonomy=config['taxonomy']),
        diamond = expand("{diamond}", diamond=config['diamond']),
    conda:
        "conda_envs/CAT.yaml"
    threads:
        24
    resources:
        mem_gb=256,
        time='24:00:00'
    benchmark:
        "3_Outputs/0_Logs/{sample}_CAT.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_CAT.log"
    message:
        "Using CAT to assign taxonomy to {wildcards.sample}'s contigs"
    shell:
        """
        CAT contigs \
            -c {input.contigs} \
            -o {params.group} \
            -d {params.database} \
            -t {params.taxonomy} \
            --path_to_diamond {params.diamond} \
            -n {threads} \
            --index_chunks 1
        
        
        #Assign tax to main output
        CAT add_names \
            -i {output.classif} \
            -t {params.taxonomy} \
            --only_official \
            -o {output.final_output}


        #Compress
        pigz -p {threads} {output.faa}
        pigz -p {threads} {output.gff}
        """