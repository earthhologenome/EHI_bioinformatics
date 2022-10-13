################################################################################
################################################################################
################################################################################
# Snakefile for individual assembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

configfile: "0_Code/2_Assembly_Binning_config.yaml"

### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_M_1.fastq.gz", "")
            for fn in glob(f"2_Reads/4_Host_removed/*_1.fastq.gz")]

print("Detected the following samples:")
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
        "3_Outputs/assembly_summary.tsv",
        "3_Outputs/5_Refined_Bins/All_bins.stats"

################################################################################
### Perform assembly on each sample
rule Assembly:
    input:
        r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq.gz",
        r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq.gz",
    output:
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
    params:
        workdir = "3_Outputs/2_Assemblies/{sample}",
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        48
    resources:
        mem_gb=256
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
            reformat.sh \
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
        "2_Assembly_Binning.yaml"
    threads:
        10
    resources:
        mem_gb=80
    message:
        "Running QUAST on {wildcards.sample} assembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.assembly}

        # Parse select metrics for final report
        grep N50 {output.report}/report.tsv | cut -f2 > {output.report}/n50.tsv
        grep L50 {output.report}/report.tsv | cut -f2 > {output.report}/l50.tsv
        grep "# contigs" {output.report}/report.tsv | cut -f2 > {output.report}/ncontigs.tsv
        grep "Largest contig" {output.report}/report.tsv | cut -f2 > {output.report}/largestcontig.tsv
        grep "Total length" {output.report}/report.tsv | cut -f2 > {output.report}/totallength.tsv

        # paste into a single table
        paste {output.report}/n50.tsv \
            {output.report}/l50.tsv \
            {output.report}/ncontigs.tsv \
            {output.report}/largestcontig.tsv \
            {output.report}/totallength.tsv > {output.report}/{wildcards.sample}_assembly_report.tsv
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
        "2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_indexing.log"
    message:
        "Indexing {wildcards.sample} assembly using Bowtie2"
    shell:
        """
        # Index the assembly
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
        r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq.gz",
        r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=80
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
### Bin each sample's contigs using MetaBAT2
rule metabat2:
    input:
        bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    output:
        metabat2_depths = "3_Outputs/4_Binning/{sample}/{sample}_metabat_depth.txt",
        metabat2 = directory("3_Outputs/4_Binning/{sample}/metabat2_bins")
    params:
        minlength = expand("{minlength}", minlength=config['minlength'])
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        48
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/{sample}_metabat2_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_metabat2_binning.log"
    message:
        "Binning {wildcards.sample} contigs with metabat2"
    shell:
        """
        # Create contig depth file
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.metabat2_depths} {input.bam}

        # Run metabat2
        metabat2 \
            -i {input.assembly} \
            -a {output.metabat2_depths} \
            -o {output.metabat2}/metabat2_bin \
            -m {params.minlength} \
            -t {threads} --unbinned
        """
################################################################################
### Bin each sample's contigs using SemiBin
rule semibin:
    input:
        bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
    output:
        semibin = directory("3_Outputs/4_Binning/{sample}/semibin_bins")
    params:
        env = expand("{env}", env=config['env'])
    conda:
        "semibin.yaml"
    threads:
        48
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/{sample}_semibin_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_semibin_binning.log"
    message:
        "Binning {wildcards.sample} contigs with semibin"
    shell:
        """
        # Run semibin
        SemiBin single_easy_bin \
            --environment {params.env} \
            -i {input.assembly} \
            -o {output.semibin} \
            -b {input.bam} \
            -t {threads}
        """
################################################################################
### Bin each sample's contigs using rosella
rule rosella:
    input:
        bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam",
        assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta",
        metabat2_depths = "3_Outputs/4_Binning/{sample}/{sample}_metabat_depth.txt"
    output:
        rosella = directory("3_Outputs/4_Binning/{sample}/rosella")
    params:
        minlength = expand("{minlength}", minlength=config['minlength'])
    conda:
        "rosella.yaml"
    threads:
        48
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/{sample}_rosella_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_rosella_binning.log"
    message:
        "Binning {wildcards.sample} contigs with rosella"
    shell:
        """
        # Run rosella
        rosella bin \
            -r {input.assembly} \
            --coverage-values {input.metabat2_depths} \
            -o {output.rosella} \
            -t {threads}
        """
################################################################################
### Bin each sample's contigs using binny
## N.B, not sure if it's possible to call binny like other binners...
## I've opened an issue here: https://github.com/a-h-b/binny/issues
# rule binny:
#     input:
#         bam = "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam",
#         assembly = "3_Outputs/2_Assemblies/{sample}_contigs.fasta"
#     output:
#         binny = directory("3_Outputs/4_Binning/{sample}/binny_bins")
#     params:
#         minlength = expand("{minlength}", minlength=config['minlength'])
#     conda:
#         "binny.yaml"
#     threads:
#         48
#     resources:
#         mem_gb=180
#     benchmark:
#         "3_Outputs/0_Logs/{sample}_metabat2_binning.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{sample}_metabat2_binning.log"
#     message:
#         "Binning {wildcards.sample} contigs with metabat2"
#     shell:
#         """
#         # Run metabat2
#         metabat2 \
#             -i {input.assembly} \
#             -o {output.binny} \
#             -m {params.minlength} \
#             -t {threads} --unbinned
#         """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        a = "3_Outputs/4_Binning/{sample}/rosella_bins",
        b = "3_Outputs/4_Binning/{sample}/semibin_bins",
        c = "3_Outputs/4_Binning/{sample}/metabat2_bins",
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
        "2_MetaWRAP.yaml"
    threads:
        48
    resources:
        mem_gb=256
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_bin_refinement.log"
    message:
        "Refining {wildcards.sample} bins with MetaWRAP's bin refinement module"
    shell:
        """
        # Setup checkM path
        export checkmdb={config[checkmdb]}
        printf $checkmdb | checkm data setRoot

        metawrap bin_refinement \
            -m {params.memory} \
            -t {threads} \
            -o {params.outdir} \
            -A {input.a} \
            -B {input.b} \
            -C {input.c} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        cp {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        cp {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/{params.sample}_bin/g' {output.stats}
        sed -i'' 's/bin/{params.sample}_bin/g' {output.contigmap}
        for bin in {params.bindir}/*.fa;
            do mv $bin ${{bin/bin./{params.sample}_bin.}};
                done

        # Compress, clean outputs:
        rm -r {params.binning_wfs}
#        rm -r {params.refinement_wfs}
        rm {input.a}/*.fa
        rm {input.b}/*.fa
        rm {input.c}/*.fa

        pigz -p {threads} {params.outdir}/*_bins/*.fa
        """
################################################################################
### Combine metawrap stats for dereplication
rule reformat_metawrap:
    input:
        expand("3_Outputs/5_Refined_Bins/{sample}_metawrap_70_10_bins.stats", sample=SAMPLE)
    output:
        stats = "3_Outputs/5_Refined_Bins/All_bins.stats",
    params:
        all_folder = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins",
        wd = "3_Outputs/5_Refined_Bins",
        stats_no_header = "3_Outputs/5_Refined_Bins/All_bins_no_header.stats",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        1
    resources:
        mem_gb=24
    message:
        "Reformatting metaWRAP outputs"
    shell:
        """
        # Copy each sample's bins to a single folder
        mkdir -p {params.all_folder}
        cp {params.wd}/*/metawrap_70_10_bins/* {params.all_folder}

        # Setup headers for combined metawrap file:
        echo -e genome'\t'completeness'\t'contamination'\t'GC'\t'lineage'\t'N50'\t'size'\t'binner > {params.wd}/header.txt

        #Cat the bin info from each group together
        for i in {params.wd}/*.stats;
            do grep -v 'contamination' $i >> {params.stats_no_header};
                done
        cat {params.wd}/header.txt {params.stats_no_header} > {output.stats}

        #Format for dRep input
        cut -f1,2,3 --output-delimiter=, {output.stats} | sed 's/,/.fa,/' | sed 's/genome.fa/bin/' > {params.wd}/All_bins_dRep.csv

        #Print the number of MAGs to a file for combining with the assembly report
        for sample in {params.wd}/*;
            do ls -l $sample/metawrap_70_10_bins/*.fa.gz | wc -l > $sample_bins.tsv;
        done

        # Clean up
        rm {params.stats_no_header}
        rm {params.wd}/header.txt
        """
################################################################################
### Calculate the number of reads that mapped to assemblies
rule coverM_assembly:
    input:
        "3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam"
    output:
        "3_Outputs/6_CoverM/{sample}_coverM_rel_abun.txt"
    params:
        sample = "{sample}"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=45
    log:
        "3_Outputs/0_Logs/{sample}_coverM_assembly.log"
    message:
        "Calculating assembly mapping rate for {wildcards.sample} with CoverM"
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
### Generate output summary table
rule generate_summary:
    input:
        expand("3_Outputs/6_CoverM/{sample}_coverM_rel_abun.txt", sample=SAMPLE),
        "3_Outputs/5_Refined_Bins/All_bins.stats"
    output:
        "3_Outputs/assembly_summary.tsv"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        1
    resources:
        mem_gb=16
    log:
        "3_Outputs/0_Logs/summarise_assembly.log"
    message:
        "Creating final summary table"
    shell:
        """
        Create the final output summary table
        #parse QUAST outputs for assembly stats
        echo -e "N50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\taseembly_mapping_percent" > headers.tsv
        cat 2_Assemblies/*_QUAST/*_assembly_report.tsv > temp_report.tsv


        #Create sampleid column
        for sample in 2_Assemblies/*_QUAST;
            do echo ${{sample/_QUAST/}} >> sampleids.tsv;
        done

        paste sampleids.tsv temp_report.tsv > temp2_report.tsv

        #Add in the # of bins
        cat *_bins.tsv > number_bins.tsv
        paste temp2_report.tsv number_bins.tsv > temp3_report.tsv

        #Add in the % mapping to assembly stats
        for sample in 3_Outputs/6_CoverM/*_coverM_rel_abun.txt;
            do sed -n 3p $sample | cut -f2 > $sample_relabun.tsv;
        done

        cat *_relabun.tsv > all_relabun.tsv

        paste temp3_report.tsv all_relabun.tsv > temp4_report.tsv

        #Combine them into the final assembly report
        cat headers.tsv temp4_report.tsv > {output}

        #Clean up
        rm headers.tsv && rm temp_report.tsv && rm temp2_report.tsv && rm *_bins.tsv
        rm *_relabun.tsv
        """