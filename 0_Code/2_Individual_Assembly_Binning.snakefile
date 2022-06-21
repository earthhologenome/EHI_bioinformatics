################################################################################
################################################################################
################################################################################
# Snakefile for individual assembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

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
        "3_Outputs/6_CoverM/coverM_assemblies_rel_abun.txt",
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
        40
    resources:
        mem_gb=180
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
        "2_Assembly_Binning.yaml"
    threads:
        20
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
        40
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
        r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq.gz",
        r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        10
    resources:
        mem_gb=45
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
        "2_MetaWRAP.yaml"
    threads:
        40
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/{sample}_assembly_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_assembly_binning.log"
    message:
        "Binning {wildcards.sample} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fastq/assembly files to trick metaWRAP into running without mapping
        mkdir -p {params.outdir}/work_files

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
            -l 1500 \
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
        "2_MetaWRAP.yaml"
    threads:
        40
    resources:
        mem_gb=180
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
        sed -i'' '2,$s/bin/{params.sample}_bin/g' {output.stats}
        sed -i'' 's/bin/{params.sample}_bin/g' {output.contigmap}
        for bin in {params.bindir}/*.fa;
            do mv $bin ${{bin/bin./{params.sample}_bin.}};
                done

        # Compress, clean outputs:
        rm -r {params.binning_wfs}
        rm -r {params.refinement_wfs}
        rm {input.concoct}/*.fa
        rm {input.maxbin2}/*.fa
        rm {input.metabat2}/*.fa

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
        sample = lambda wildcards: {sample}
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

        # Clean up
        rm {params.stats_no_header}
        rm {params.wd}/header.txt
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        expand("3_Outputs/3_Assembly_Mapping/BAMs/{sample}.bam", sample=SAMPLE)
    output:
        "3_Outputs/6_CoverM/coverM_assemblies_rel_abun.txt"
    params:
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=45
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
