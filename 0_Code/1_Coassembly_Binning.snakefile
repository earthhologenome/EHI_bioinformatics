################################################################################
################################################################################
################################################################################
### Snakefile for EHI's coassembly and binning pipeline                      ###
### Raphael Eisenhofer 01/2022                                               ###
###                                                                          ###
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

GROUP = [os.path.basename(dir)
         for dir in glob(f"2_Reads/3_Host_removed/*")]

SAMPLE = [os.path.relpath(fn, "2_Reads/3_Host_removed").replace("_non_host_1.fastq.gz", "")
            for group in GROUP
            for fn in glob(f"2_Reads/3_Host_removed/{group}/*_1.fastq.gz")]

print("Detected these sample groups:")
print(GROUP)
print("Detected the following samples:")
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("2_Reads/1_Trimmed/{sample}_trimmed_1.fastq.gz", sample=SAMPLE),
        expand("2_Reads/1_Trimmed/{sample}_trimmed_2.fastq.gz", sample=SAMPLE),
        expand("3_Outputs/1_QC/1_Host_BAMs/{sample}_host.bam", sample=SAMPLE),
        "3_Outputs/1_QC/2_CoverM/coverM_mapped_host.tsv",
        expand("3_Outputs/3_Coassembly_Mapping/BAMs/{group}_coverM.txt", group=GROUP)
################################################################################
### Pull samples from ENA            *** IN PROGRESS ***
### Probably use grabseqs: https://github.com/louiejtaylor/grabseqs
# rule pull_samples:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
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
### Download host reference genome        *** IN PROGRESS ***
###
# rule pull_reference:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:

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
### Create figures/table from host CoverM output        *** IN PROGRESS ***
###
# rule create_host_mapping_reports:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Perform Coassemblies on each sample group
rule Coassembly:
    input:
        reads = "2_Reads/3_Host_removed/{group}"
    output:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        r1_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_1.fastq.gz"),
        r2_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_2.fastq.gz")
    params:
        workdir = "3_Outputs/2_Coassemblies/{group}",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly.log"
    message:
        "Coassembling {wildcards.group} using metaspades"
    shell:
        """
        # Concatenate reads from the same group for Coassembly
        cat {input.reads}/*_1.fastq.gz > {output.r1_cat}
        cat {input.reads}/*_1.fastq.gz > {output.r2_cat}

        # Run metaspades
        metaspades.py \
            -t {threads} \
            -k 21,33,55,77,99 \
            -1 {output.r1_cat} -2 {output.r2_cat} \
            -o {params.workdir}
        2> {log}

        # Move the Coassembly to final destination
        mv {params.workdir}/scaffolds.fasta {output.Coassembly}
        """
################################################################################
### Create QUAST reports of coassemblies
rule QUAST:
    input:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        report = "3_Outputs/2_Coassemblies/{group}_QUAST/report.html",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        20
    message:
        "Running QUAST on {wildcards.group} coassembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.Coassembly}
        """
################################################################################
### Map reads to the coassemblies
rule Coassembly_index:
    input:
        report = "3_Outputs/2_Coassemblies/{group}_QUAST/report.html"
    output:
        bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l",
    params:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_indexing.log"
    message:
        "Indexing {wildcards.group} coassembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.Coassembly} {params.Coassembly} \
        &> {log}
        """
################################################################################
### Map reads to the coassemblies
rule Coassembly_mapping:
    input:
        bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l"
    output:
        mapped_bam = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}"
    params:
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        read_dir = "2_Reads/3_Host_removed/{group}"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_mapping.log"
    message:
        "Mapping {wildcards.group} samples to coassembly using Bowtie2"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        for fq1 in {params.read_dir}/*_1.fastq.gz; do \
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.assembly} \
            -1 $fq1 \
            -2 ${{fq1/_1.fastq.gz/_2.fastq.gz}} \
        | samtools view -@ {threads} -o {output.mapped_bam}/${{fq1/_1.fastq.gz/.bam}} -; done
        """
################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Coassembly_Mapping/BAMs/{group}"
    output:
        concoct = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/metabat2_bins",
    params:
        outdir = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/Binning",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        memory = "16"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_binning.log"
    message:
        "Binning {wildcards.group} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input}/*.bam; do ln -s $bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            -a {params.assembly} \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq
        """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        concoct = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/3_Coassembly_Mapping/Binning/{group}/metabat2_bins",
    output:
        stats = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.contigs"
    params:
        outdir = "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}",
        memory = "16",
        group = "{group}"
    conda:
        "2_MetaWRAP.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Refining {wildcards.group} bins with MetaWRAP's bin refinement module"
    shell:
        """
        metawrap bin_refinement \
            -m {params.memory} \
            -t {threads} \
            -o {params.outdir} \
            -A {input.concoct} \
            -B {input.maxbin2} \
            -C {input.metabat2} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/bin_{params.group}/g' {output.stats}
        sed -i'' 's/bin/bin_{params.group}/g' {output.contigmap}
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        "3_Outputs/3_Coassembly_Mapping/Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats"
    output:
        "3_Outputs/3_Coassembly_Mapping/BAMs/{group}_coverM.txt"
    params:
        mapped_bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        memory = "16",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Calculating coassembly mapping rate for {wildcards.group} with CoverM"
    shell:
        """
        coverm genome \
            -b {params.mapped_bams}/*.bam \
            --genome-fasta-files {params.assembly} \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
################################################################################
### Create figures/table from MAG mapping results        *** IN PROGRESS ***
###
# rule create_MAG_mapping_reports:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Create figures/table QUAST outputs        *** IN PROGRESS ***
###
# rule create_QUAST_reports:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Upload outputs to ENA        *** IN PROGRESS ***
###
# rule upload_to_ENA:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################

onsuccess:
    shell("""
            mail -s "workflow completed" raph.eisenhofer@gmail.com < {log}

            #Clean up files
            rm 3_Outputs/1_QC/1_BAMs/*/*.bam
          """)
