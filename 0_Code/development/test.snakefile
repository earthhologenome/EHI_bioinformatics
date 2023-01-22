import os
from glob import glob

GROUPS = {}

# The directory where the 'groupN' folders are located
parent_dir = '2_Reads/4_Host_removed'

# Iterate through all the directories in the specified directory
for group_name in os.listdir(parent_dir):
    # Construct the full path to the current group directory
    group_path = os.path.join(parent_dir, group_name)
    # Check if the current item is a directory
    if os.path.isdir(group_path):
        # Initialize an empty list to store file names
        samples = []
        # Iterate through all the files in the current directory
        for sample_name in os.listdir(group_path):
            # Construct the full path to the current sample file
            sample_path = os.path.join(group_path, sample_name)
            # Check if the current item is a file
            if os.path.isfile(sample_path) and sample_name.endswith("_M_1.fastq.gz"):
                # Append the file name to the list of samples
                samples.append(os.path.basename(sample_path).replace("_M_1.fastq.gz", ""))
        # Add the group name and list of samples to the dictionary
        GROUPS[group_name] = samples



print(GROUPS)

configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"
  

rule all:
    input:
       expand("3_Outputs/{group}_coassembly_summary.tsv", group=GROUPS.keys()),



rule Coassembly:
    input:
        reads = "2_Reads/4_Host_removed/{group}"
    output:
        coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    params:
        workdir = "3_Outputs/2_Coassemblies/{group}",
        r1_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_1.fastq.gz"),
        r2_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_2.fastq.gz"),
        assembler = expand("{assembler}", assembler=config['assembler']),
    shell:
        """
        # Set up assembler variable from config file
        export assembler={config[assembler]}

        if [ "$assembler" == "metaspades" ]
        then

        # Concatenate reads from the same group for Coassembly
        cat {input.reads}/*_1.fastq.gz > {params.r1_cat}
        cat {input.reads}/*_2.fastq.gz > {params.r2_cat}

        # Run metaspades
            metaspades.py \
                -t {threads} \
                -k 21,33,55,77,99 \
                -1 {params.r1_cat} -2 {params.r2_cat} \
                -o {params.workdir}
                2> {log}

        # Remove contigs shorter than 1,500 bp
            reformat.sh \
                in={params.workdir}/scaffolds.fasta \
                out={output.coassembly} \
                minlength=1500

        else

        # Set up input reads variable for megahit
        R1=$(for i in 2_Reads/4_Host_removed/{wildcards.group}/*_1.fastq.gz; do echo $i | tr '\n' ,; done)
        R2=$(for i in 2_Reads/4_Host_removed/{wildcards.group}/*_2.fastq.gz; do echo $i | tr '\n' ,; done)

        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 $R1 -2 $R2 \
                -f \
                -o {params.workdir}
                2> {log}

        # Move the Coassembly to final destination
            mv {params.workdir}/final.contigs.fa {output.coassembly}

        # Reformat headers
            sed -i 's/ /-/g' {output.coassembly}

        fi        
        """

rule mapping:
    input:
        contigs = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        "3_Outputs/3_BAMs/{group}/{sample}.bam"  
    params:
        r1 = "2_Reads/4_Host_removed/{group}/{sample}_M_1.fastq.gz",
        r2 = "2_Reads/4_Host_removed/{group}/{sample}_M_2.fastq.gz",
        threads = 24
    shell:
            """
            # Map reads to catted reference using Bowtie2
            bowtie2 \
                --time \
                --threads {threads} \
                -x {input.contigs} \
                -1 {params.r1} \
                -2 {params.r2} \
            | samtools sort -@ {threads} -o {output}
            """

rule summary:
    input:
       bams = lambda wildcards: ["3_Outputs/3_BAMs/{}/{}.bam/".format(wildcards.group, sample) for sample in GROUPS[wildcards.group]],
    output:
        "3_Outputs/{group}_coassembly_summary.tsv"
    shell:
        "create_summary.py {input.bams} > {output}"