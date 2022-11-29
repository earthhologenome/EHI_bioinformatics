################################################################################
################################################################################
################################################################################
# Snakefile for individual assembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"

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
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=128,
        time='12:00:00'
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
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=64,
        time='00:15:00'
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
        grep "# contigs (>= 0 bp)" {output.report}/report.tsv | cut -f2 > {output.report}/ncontigs.tsv
        grep "Largest contig" {output.report}/report.tsv | cut -f2 > {output.report}/largestcontig.tsv
        grep "Total length (>= 0 bp)" {output.report}/report.tsv | cut -f2 > {output.report}/totallength.tsv

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
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=128,
        time='02:00:00'
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
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=128,
        time='04:00:00'
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
        "conda_envs/2_MetaWRAP.yaml"
    threads:
        24
    resources:
        mem_gb=128,
        time='08:00:00'
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
        ln -sf `pwd`/{input} {params.outdir}/work_files/$(basename {input})

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
        stats = "3_Outputs/5_Refined_Bins/{sample}/{sample}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/5_Refined_Bins/{sample}/{sample}_metawrap_70_10_bins.contigs",
        outdir = directory("3_Outputs/5_Refined_Bins/{sample}")
    params:
        outdir = "3_Outputs/5_Refined_Bins/{sample}",
        bindir = "3_Outputs/5_Refined_Bins/{sample}/metawrap_70_10_bins",
        binning_wfs = "3_Outputs/4_Binning/{sample}/work_files",
        refinement_wfs = "3_Outputs/5_Refined_Bins/{sample}/work_files",
        memory = "180",
        sample = "{sample}"
    conda:
        "conda_envs/2_MetaWRAP.yaml"
    threads:
        24
    resources:
        mem_gb=256,
        time='06:00:00'
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
            -A {input.concoct} \
            -B {input.maxbin2} \
            -C {input.metabat2} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
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
#        expand("3_Outputs/5_Refined_Bins/{sample}/{sample}_metawrap_70_10_bins.stats", sample=SAMPLE)
        expand("3_Outputs/5_Refined_Bins/{sample}", sample=SAMPLE)
    output:
        stats = "3_Outputs/5_Refined_Bins/All_bins.stats",
    params:
        all_folder = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins",
        wd = "3_Outputs/5_Refined_Bins",
        stats_no_header = "3_Outputs/5_Refined_Bins/All_bins_no_header.stats",
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        1
    resources:
        mem_gb=24,
        time='00:05:00'
    message:
        "Reformatting metaWRAP outputs"
    shell:
        """
        #Print the number of MAGs to a file for combining with the assembly report
        for sample in {input};
            do ls -l $sample/metawrap_70_10_bins/*.fa.gz | wc -l > $(basename $sample)_bins.tsv;
        done

        # Copy each sample's bins to a single folder
        mkdir -p {params.all_folder}
        for i in {input};
            do cp $i/metawrap_70_10_bins/*.gz {params.all_folder}
        done


        # Setup headers for combined metawrap file:
        echo -e genome'\t'completeness'\t'contamination'\t'GC'\t'lineage'\t'N50'\t'size'\t'binner > {params.wd}/header.txt

        #Cat the bin info from each group together
         for i in {input};
            do grep -v 'contamination' $i/*_metawrap_70_10_bins.stats >> {params.stats_no_header};
                done
        cat {params.wd}/header.txt {params.stats_no_header} > {output.stats}

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
        "3_Outputs/6_Assembly_CoverM/{sample}_coverM_rel_abun.txt"
    params:
        sample = "{sample}"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=45,
        time='02:00:00'
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
        expand("3_Outputs/6_Assembly_CoverM/{sample}_coverM_rel_abun.txt", sample=SAMPLE),
        "3_Outputs/5_Refined_Bins/All_bins.stats"
    output:
        "3_Outputs/assembly_summary.tsv"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        1
    resources:
        mem_gb=16,
        time='00:05:00'
    log:
        "3_Outputs/0_Logs/summarise_assembly.log"
    message:
        "Creating final summary table"
    shell:
        """
        #Create the final output summary table
        #parse QUAST outputs for assembly stats
        echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\taseembly_mapping_percent" > headers.tsv
        cat 3_Outputs/2_Assemblies/*_QUAST/*_assembly_report.tsv > temp_report.tsv


        #Create sampleid column
        for sample in 3_Outputs/2_Assemblies/*_QUAST;
            do echo $(basename ${{sample/_QUAST/}}) >> sampleids.tsv;
        done

        paste sampleids.tsv temp_report.tsv > temp2_report.tsv

        #Add in the # of bins
        cat *_bins.tsv > number_bins.tsv
        paste temp2_report.tsv number_bins.tsv > temp3_report.tsv

        #Add in the % mapping to assembly stats
        for sample in 3_Outputs/6_Assembly_CoverM/*_coverM_rel_abun.txt;
            do sed -n 3p $sample | cut -f2 > $(basename ${{sample/_coverM_rel_abun.txt/}})_relabun.tsv;
        done

        cat *_relabun.tsv > all_relabun.tsv

        paste temp3_report.tsv all_relabun.tsv > temp4_report.tsv

        #Combine them into the final assembly report
        cat headers.tsv temp4_report.tsv > {output}

        #Clean up
        rm *.tsv
        """