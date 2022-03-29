################################################################################
################################################################################
################################################################################
# Snakefile for quality controlling (trimming/mapping to host) reads and assembly
# Raphael Eisenhofer 10/2021
#
################################################################################
################################################################################
################################################################################

### Setup group and sample wildcards:
import os
from glob import glob

# GROUP = [os.path.basename(dir)
#          for dir in glob(f"2_Reads/0_Untrimmed/*")]
# SAMPLE = [os.path.relpath(fn, "2_Reads/0_Untrimmed/").replace("_1.fastq.gz", "")
#             for group in GROUP
#             for fn in glob(f"2_Reads/0_Untrimmed/{group}/*_1.fastq.gz")]

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"2_Reads/0_Untrimmed/*_1.fastq.gz")]

# print(GROUP)
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz", sample=SAMPLE),
        expand("2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz", sample=SAMPLE),
        expand("3_Outputs/1_QC/1_Host_BAMs/{sample}_host.bam", sample=SAMPLE),
        "3_Outputs/1_QC/2_CoverM/coverM_mapped_host.tsv"
################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i = "2_Reads/0_Untrimmed/{sample}_1.fastq.gz",
        r2i = "2_Reads/0_Untrimmed/{sample}_2.fastq.gz"
    output:
        r1o = "2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2o = "2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz",
        fastp_html = "2_Reads/2_fastp_results/{sample}.html",
        fastp_json = "2_Reads/2_fastp_results/{sample}.json"
    conda:
        "1_QC.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_fastp.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --adapter_sequence_r2 CTGTCTCTTATACACATCT \
        &> {log}
        """
################################################################################
## Index host genomes:
rule index_ref:
    input:
        "1_References"
    output:
        bt2_index = "1_References/CattedRefs.fna.gz.rev.2.bt2l",
        catted_ref = "1_References/CattedRefs.fna.gz"
    conda:
        "1_QC.yaml"
    threads:
        40
    log:
        "3_Outputs/0_Logs/host_genome_indexing.log"
    message:
        "Concatenating and indexing host genomes with Bowtie2"
    shell:
        """
        # Concatenate input reference genomes
        cat {input}/*.gz > {input}/CattedRefs.fna.gz

        # Index catted genomes
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.catted_ref} {output.catted_ref} \
        &> {log}
        """
################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i = "2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2i = "2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz",
        catted_ref = "1_References/CattedRefs.fna.gz",
        bt2_index = "1_References/CattedRefs.fna.gz.rev.2.bt2l"
    output:
        all_bam = "3_Outputs/1_QC/1_BAMs/{sample}.bam",
        host_bam = "3_Outputs/1_QC/1_Host_BAMs/{sample}_host.bam",
        non_host_r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        non_host_r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz",
    conda:
        "1_QC.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_mapping.log"
    message:
        "Mapping {wildcards.sample} reads to host genomes"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.catted_ref} \
            -1 {input.r1i} \
            -2 {input.r2i} \
        | samtools view -b -@ {threads} - | samtools sort -@ {threads} -o {output.all_bam} - &&

        # Split extract non-host reads
        samtools view -b -f12 -@ {threads} {output.all_bam} \
        | samtools fastq -@ {threads} -c 6 -1 {output.non_host_r1} -2 {output.non_host_r2} - &&

        # Send host reads to BAM
        samtools view -b -F12 -@ {threads} {output.all_bam} \
        | samtools sort -@ {threads} -o {output.host_bam} -
        """
################################################################################
### Calculate % of each sample's reads mapping to host genome/s
rule coverM:
    input:
        expand("3_Outputs/1_QC/1_BAMs/{sample}.bam", sample=SAMPLE)
    output:
        "3_Outputs/1_QC/2_CoverM/coverM_mapped_host.tsv"
    params:
        assembly = "1_References/CattedRefs.fna.gz"
    conda:
        "1_QC.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/coverM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/coverM.log"
    message:
        "Calculating percentage of reads mapped to host genome/s using coverM"
    shell:
        """
        #Calculate % mapping to host using coverM
        coverm genome \
            -b {input} \
            -s _ \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
################################################################################
onsuccess:
    shell("""
            mail -s "workflow completed" raph.eisenhofer@gmail.com < {log}

            #Clean up files
            rm 3_Outputs/1_QC/1_BAMs/*/*.bam
          """)
