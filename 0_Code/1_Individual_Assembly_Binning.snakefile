################################################################################
################################################################################
################################################################################
### Snakefile for EHI's individual assembly and binnsing pipeline            ###
### Raphael Eisenhofer 01/2022                                               ###
###                                                                          ###
################################################################################
################################################################################
################################################################################

### Setup group and sample wildcards:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"2_Reads/0_Untrimmed/*_1.fastq.gz")]

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
        "3_Outputs/6_CoverM/coverM_assemblies_rel_abun.txt",
        "3_Outputs/5_Refined_Bins/All_bins.stats"
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
### Create figures/table from CoverM output        *** IN PROGRESS ***
###
# rule create_host_mapping_reports:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Perform assembly on each sample
rule Assembly:
    input:
        r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz",
    output:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
    params:
        workdir = "3_Outputs/2_Assemblies/{sample}",
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly.log"
    message:
        "Assembling {wildcards.sample} using {params.assembler}"
    shell:
        """
        # Set up assembler variable from config file
        export assembler={config[assembler]}

        if [ "$assembler" == "metaspades" ]
        then
        # Run metaspades
            metaspades.py \
                -t {threads} \
                -k 21,33,55,77,99 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.workdir}
                2> {log}

        # Remove contigs shorter than 1,500 bp
            reformat.sh
                in={params.workdir}/scaffolds.fasta \
                out={output.assembly} \
                minlength=1500

        else
        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -o {params.workdir}
                2> {log}

        # Move the Coassembly to final destination
            mv {params.workdir}/final.contigs.fa {output.assembly}

        # Reformat headers
            sed -i 's/ /-/g' {output.assembly}

        fi
        """
################################################################################
### Create QUAST reports of assemblies
rule QUAST:
    input:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    output:
        report = directory("3_Outputs/2_Assemblies/{sample}_QUAST")
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        8
    message:
        "Running QUAST on {wildcards.sample} assembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.assembly}
        """
################################################################################
### Index each sample's assembly
rule assembly_index:
    input:
        report = "3_Outputs/2_Assemblies/{sample}_QUAST"
    output:
        bt2_index = "3_Outputs/2_Assemblies/{sample}_contigs.fasta.rev.2.bt2l"
    params:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.log"
    message:
        "Indexing {wildcards.sample} assembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.assembly} {params.assembly} \
        &> {log}
        """
################################################################################
### Map a sample's reads to it corresponding assembly
rule assembly_mapping:
    input:
        bt2_index = "3_Outputs/2_Assemblies/{sample}_contigs.fasta.rev.2.bt2l"
    output:
        mapped_bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    params:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        r1 = "2_Reads/3_Host_removed/{sample}_non_host_1.fastq.gz",
        r2 = "2_Reads/3_Host_removed/{sample}_non_host_2.fastq.gz"
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_mapping.log"
    message:
        "Mapping {wildcards.sample} to its assembly using Bowtie2"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.assembly} \
            -1 {params.r1} \
            -2 {params.r2} \
        | samtools sort -@ {threads} -o {output.mapped_bam}
        """
################################################################################
### Bin each sample's contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        concoct = directory("3_Outputs/4_Binning/{sample}/concoct_bins"),
        maxbin2 = directory("3_Outputs/4_Binning/{sample}/maxbin2_bins"),
        metabat2 = directory("3_Outputs/4_Binning/{sample}/metabat2_bins")
    params:
        outdir = "3_Outputs/4_Binning/{sample}",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        basename = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}",
        memory = "180"
    conda:
        "0_EHI_Bioinformatics_MetaWRAP_env.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_binning.log"
    message:
        "Binning {wildcards.sample} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_1.fastq)
        echo "@" > {params.outdir}/work_files/$(basename {params.basename}_2.fastq)

        #Symlink BAMs for metaWRAP
        ln -s `pwd`/{input} {params.outdir}/work_files/$(basename {input})

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
        concoct = "3_Outputs/4_Binning/{sample}/concoct_bins",
        maxbin2 = "3_Outputs/4_Binning/{sample}/maxbin2_bins",
        metabat2 = "3_Outputs/4_Binning/{sample}/metabat2_bins",
    output:
        stats = "3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.contigs"
    params:
        outdir = "3_Outputs/5_Refined_Bins/{sample}",
        bindir = "3_Outputs/5_Refined_Bins/{sample}/metawrap_70_10_bins",
        binning_wfs = "3_Outputs/4_Binning/{sample}/work_files",
        refinement_wfs = "3_Outputs/5_Refined_Bins/{sample}/work_files",
        memory = "180",
        sample = "{sample}"
    conda:
        "0_EHI_Bioinformatics_MetaWRAP_env.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.log"
    message:
        "Refining {wildcards.sample} bins with MetaWRAP's bin refinement module"
    shell:
        """
        # Setup checkM path
        printf "/home/projects/ku-cbd/people/rapeis/0_DBs/CHECKM" | checkm data setRoot

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
        cp {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        cp {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/bin_{params.sample}/g' {output.stats}
        sed -i'' 's/bin/bin_{params.sample}/g' {output.contigmap}
        for bin in {params.bindir}/*.fa;
            do mv $bin ${{bin/bin/bin_{params.sample}}};
                done

        # Compress, clean outputs:
        rm -r {params.binning_wfs}/*_out
        rm {params.binning_wfs}/assembly*
        rm -r {params.refinement_wfs}
        rm -r {params.outdir}/concoct_bins
        rm -r {params.outdir}/maxbin2_bins
        rm -r {params.outdir}/metabat2_bins

        pigz -p {threads} {params.bindir}/*
        """
################################################################################
### Reformat metawrap outputs
rule reformat_metawrap:
    input:
        expand("3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.stats", sample=SAMPLE)
    output:
        stats = "3_Outputs/5_Refined_Bins/All_bins.stats",
    params:
        all_folder = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins",
        wd = "3_Outputs/5_Refined_Bins",
        stats_no_header = "3_Outputs/5_Refined_Bins/All_bins_no_header.stats",
        sample = "{sample}"
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        8
    message:
        "Reformatting metaWRAP outputs"
    shell:
        """
        # Copy each sample's bins to a single folder
        mkdir -p {params.all_folder}
        cp {params.wd}/*/metawrap_70_10_bins/* {params.all_folder}

        # Setup headers for combined metawrap file:
        echo -e bin' \t 'completeness' \t 'contamination' \t 'GC' \t 'lineage' \t 'N50' \t 'size' \t 'binner > {params.wd}/header.txt

        #Cat the bin info from each group together
        grep -v 'contamination' {params.wd}/{params.sample}_metawrap_70_10_bins.stats >> {params.stats_no_header}
        cat {params.wd}/header.txt {params.stats_no_header} > {output.stats}

        # Clean up
        rm {params.stats_no_header}
        rm {params.wd}/header.txt
        """
################################################################################
### Functionally annotate MAGs        *** IN PROGRESS ***
###
# rule annotate_function:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Taxonomically annotate MAGs        *** IN PROGRESS ***
###
# rule annotate_taxonomy:
#     input:
#     output:
#     params:
#     message:
#     threads:
#     shell:
################################################################################
### Calculate the number of reads that mapped to assemblies
rule coverM_assembly:
    input:
        expand("3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam", sample=SAMPLE)
    output:
        "3_Outputs/6_CoverM/coverM_assemblies_rel_abun.txt"
    params:
    conda:
        "0_EHI_Bioinformatics_conda_env.yaml"
    threads:
        8
    benchmark:
        "3_Outputs/0_Logs/coverM_assembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/coverM_assembly.log"
    message:
        "Calculating assembly mapping rate with CoverM"
    shell:
        """
        coverm genome \
            -b {input} \
            -m relative_abundance \
            -t {threads} \
            -s _ \
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
