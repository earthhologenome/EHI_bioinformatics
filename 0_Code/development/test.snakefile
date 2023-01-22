import os
from glob import glob

GROUP = [ dir for dir in os.listdir('2_Reads/4_Host_removed/')
         if os.path.isdir(os.path.join('2_Reads/4_Host_removed/', dir)) ]

SAMPLE = [os.path.basename(fn).replace("_M_1.fastq.gz", "")
            for fn in glob(f"2_Reads/4_Host_removed/*/*_M_1.fastq.gz")]


print(GROUP)
print(SAMPLE)

configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"
  


rule all:
    input:
        expand("3_Outputs/4_Summary/{group}_coassembly_summary.tsv", group=GROUP)


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
        r1 = lambda wildcards: ["2_Reads/4_Host_removed/{group}/{sample}_M_1.fastq.gz".format(group=group, sample=sample)
                               for group in GROUP
                               for sample in SAMPLE
                               if os.path.isfile("2_Reads/4_Host_removed/{group}/{sample}_M_1.fastq.gz".format(group=group, sample=sample))],
        contigs = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        "3_Outputs/3_BAMs/{group}/{sample}.bam"   
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.contigs} \
            -1 {input.r1} \
        | samtools sort -@ {threads} -o {output}
        """

rule summary:
    input:
        expand("3_Outputs/3_BAMs/{group}/{sample}.bam", group=GROUP, sample=SAMPLE)
    output:
        expand("3_Outputs/4_Summary/{group}_coassembly_summary.tsv", group=GROUP)
    shell:
        "create_summary.py {input} > {output}"
