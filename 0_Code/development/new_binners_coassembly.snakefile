################################################################################
################################################################################
################################################################################
# Snakefile for coassembly, binning, and refinement of MAGs
# Raphael Eisenhofer 11/2021
#
################################################################################
################################################################################
################################################################################

configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"

### Setup sample wildcard:
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
#        expand("3_Outputs/6_CoverM/{group}_assembly_coverM.txt", group=GROUP)
        expand("3_Outputs/{group}_coassembly_summary.tsv", group=GROUP)

################################################################################
### Perform Coassemblies on each sample group
rule Coassembly:
    input:
        reads = "2_Reads/4_Host_removed/{group}"
    output:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
    params:
        workdir = "3_Outputs/2_Coassemblies/{group}",
        r1_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_1.fastq.gz"),
        r2_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_2.fastq.gz"),
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        48
    resources:
        mem_gb=512,
        time='36:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly.log"
    message:
        "Coassembling {wildcards.group} using {params.assembler}"
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
                out={output.Coassembly} \
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
            mv {params.workdir}/final.contigs.fa {output.Coassembly}

        # Reformat headers
            sed -i 's/ /-/g' {output.Coassembly}

        fi
        """
################################################################################
### Create QUAST reports of coassemblies
rule QUAST:
    input:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    output:
        report = directory("3_Outputs/2_Coassemblies/{group}_QUAST"),
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        8
    resources:
        mem_gb=90,
        time='01:00:00'
    message:
        "Running -QUAST on {wildcards.group} coassembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output.report} \
            --threads {threads} \
            {input.Coassembly}

        # # Rename QUAST files
        # for i in {output.report}/*;
        #     do mv $i {output.report}/{wildcards.group}_$(basename $i);
        #         done
        
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
              {output.report}/totallength.tsv > {output.report}/{wildcards.group}_assembly_report.tsv
        """
################################################################################
### Map reads to the coassemblies
rule Coassembly_index:
    input:
        report = "3_Outputs/2_Coassemblies/{group}_QUAST/"
    output:
        bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l",
    params:
        Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        32
    resources:
        mem_gb=180,
        time='04:00:00'
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
        directory("3_Outputs/3_Coassembly_Mapping/BAMs/{group}/Complete")
    params:
        outdir = directory("3_Outputs/3_Coassembly_Mapping/BAMs/{group}"),
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        read_dir = "2_Reads/4_Host_removed/{group}"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        32
    resources:
        mem_gb=256,
        time='24:00:00'
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
        | samtools sort -@ {threads} -o {params.outdir}/$(basename ${{fq1/_1.fastq.gz/.bam}}); done

        #Create output file for snakemake
        mkdir -p {output}
        """
################################################################################
### Bin each sample's contigs using MetaBAT2
rule metabat2:
    input:
        bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
    output:
        metabat2_depths = "3_Outputs/4_Binning/{group}/{group}_metabat_depth.txt",
        metabat2 = directory("3_Outputs/4_Binning/{group}/metabat2_bins")
    params:
        minlength = expand("{minlength}", minlength=config['minlength'])
    conda:
        "conda_envs/metabat2.yaml"
    threads:
        32
    resources:
        mem_gb=180,
        time='24:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_metabat2_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_metabat2_binning.log"
    message:
        "Binning {wildcards.group} contigs with metabat2"
    shell:
        """
        # Create contig depth file
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.metabat2_depths} {input.bams}/*.bam

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
        bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
    output:
        semibin = directory("3_Outputs/4_Binning/{group}/semibin_bins/output_recluster_bins")
    params:
        env = expand("{env}", env=config['env']),
        output = "3_Outputs/4_Binning/{group}/semibin_bins/"
    conda:
        "conda_envs/semibin.yaml"
    threads:
        32
    resources:
        mem_gb=180,
        time='54:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_semibin_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_semibin_binning.log"
    message:
        "Binning {wildcards.group} contigs with semibin"
    shell:
        """
        # Run semibin
        SemiBin single_easy_bin \
            -r /projects/mjolnir1/people/ncl550/0_software/GTDB_SemiBin \
            --self-supervised \
            --training-type self \
            -i {input.assembly} \
            -o {params.output} \
            -b {input.bams}/*.bam \
            -t {threads}
        """
################################################################################
### Bin each sample's contigs using metabinner
rule metabinner:
    input:
        bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly_dir = "3_Outputs/2_Coassemblies/{group}/",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        metabat2_depths = "3_Outputs/4_Binning/{group}/{group}_metabat_depth.txt"
    output:
        coverage_file = "3_Outputs/4_Binning/{group}/coverage_profile.tsv"
        metabinner = directory("3_Outputs/4_Binning/{group}/metabinner_bins")
    params:
        minlength = expand("{minlength}", minlength=config['minlength']),
        dir = "3_Outputs/4_Binning/{group}"
    conda:
        "conda_envs/metabinner.yaml"
    threads:
        32
    resources:
        mem_gb=180,
        time='54:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_metabinner_binning.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_metabinner_binning.log"
    message:
        "Binning {wildcards.group} contigs with metabinner"
    shell:
        """
        # Edit metabat2 coverage profile for use with metabinner
        cat {input.metabat2_depths} | cut -f -1,4- > {output.coverage_file}

        # Generate contig kmer profiles
        cd {input.assembly_dir}
        python gen_kmer.py *contigs.fasta 1500 4

        # Run metabinner
        bash run_metabinner.sh \
            -a {input.assembly} \
            -d {output.coverage_file} \
            -k {input.assembly_dir}/*kmer_4_f1500.csv \
            -o {output.metabinner} \
            -t {threads}

        """
# ################################################################################
# ### Bin each sample's contigs using rosella
# rule rosella:
#     input:
#         bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
#         assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
#         metabat2_depths = "3_Outputs/4_Binning/{group}/{group}_metabat_depth.txt"
#     output:
#         rosella = directory("3_Outputs/4_Binning/{group}/rosella_bins")
#     params:
#         minlength = expand("{minlength}", minlength=config['minlength']),
#         dir = "3_Outputs/4_Binning/{group}"
#     conda:
#         "rosella.yaml"
#     threads:
#         32
#     resources:
#         mem_gb=180,
#         time='54:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_rosella_binning.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_rosella_binning.log"
#     message:
#         "Binning {wildcards.group} contigs with rosella"
#     shell:
#         """
#         # Run rosella
#         rosella recover \
#             -r {input.assembly} \
#             --coverage-values {input.metabat2_depths} \
#             -o {output.rosella} \
#             -t {threads}

#         # Move misc files from rosella output in case they interfere with refinement
#         mv {output.rosella}/*.png {params.dir}
#         mv {output.rosella}/*.tsv {params.dir}
#         mv {output.rosella}/*.json {params.dir}
#         """
################################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        a = "3_Outputs/4_Binning/{group}/metabinner_bins",
        b = "3_Outputs/4_Binning/{group}/semibin_bins/output_recluster_bins",
        c = "3_Outputs/4_Binning/{group}/metabat2_bins",    
    output:
        stats = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats",
        contigmap = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.contigs"
    params:
        refinement_wfs = "3_Outputs/5_Refined_Bins/{group}/work_files",
        outdir = "3_Outputs/5_Refined_Bins/{group}",
        memory = "256",
        group = "{group}"
    conda:
        "conda_envs/2_MetaWRAP.yaml"
    threads:
        32
    resources:
        mem_gb=256,
        time='36:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Refining {wildcards.group} bins with MetaWRAP's bin refinement module"
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
        mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/{params.group}_bin/g' {output.stats}
        sed -i'' 's/bin/{params.group}_bin/g' {output.contigmap}
        for bin in {params.outdir}/metawrap_70_10_bins/*.fa;
            do mv $bin ${{bin/bin./{params.group}_bin.}};
                done

        # Compress output bins
        pigz -p {threads} {params.outdir}/*bins/*.fa

        rm -r {params.refinement_wfs}
        rm {input.a}/*.fa
        rm {input.b}/*.fa
        rm {input.c}/*.fa
        """
################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats"
    output:
        "3_Outputs/6_CoverM/{group}_assembly_coverM.txt"
    params:
        mapped_bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
        assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
        binning_files = "3_Outputs/4_Binning/{group}",
        refinement_files = "3_Outputs/5_Refined_Bins/{group}",
        memory = "180",
        group = "{group}"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads:
        24
    resources:
        mem_gb=180,
        time='06:00:00'
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

        # Create directory for dereplication groups:
        mkdir -p 3_Outputs/5_Refined_Bins/dRep_groups

        #Print the number of MAGs to a file for combining with the assembly report
        ls -l {params.refinement_files}/metawrap_70_10_bins/*.fa.gz | wc -l > {params.group}_bins.tsv;
        """
################################################################################
### Generate output summary table
rule generate_summary:
    input:
        "3_Outputs/6_CoverM/{group}_assembly_coverM.txt"
    output:
        "3_Outputs/{group}_coassembly_summary.tsv"
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    params:
        group = "{group}"
    threads:
        1
    resources:
        mem_gb=16,
        time='00:05:00'
    message:
        "Creating final coassembly summary table"
    shell:
        """
        #Create the final output summary table
        #parse QUAST outputs for assembly stats
        echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\taseembly_mapping_percent" > headers.tsv
        cat 3_Outputs/2_Coassemblies/{params.group}_QUAST/{params.group}_assembly_report.tsv > {params.group}_temp_report.tsv

        # #Add in the % mapping to assembly stats
        # for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
        #     do echo $(basename ${{sample/.bam/}}) > {params.group}_$(basename "$sample")_id.tsv;
        # done

        #Add in the % mapping to assembly stats
        for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
            do echo $(basename ${{sample/.bam/}}) >> {params.group}_sample_ids.tsv;
        done

        paste {params.group}_sample_ids.tsv {params.group}_temp_report.tsv > {params.group}_temp2_report.tsv

        #Add in the # of bins
        cat {params.group}_bins.tsv > {params.group}_number_bins.tsv
        paste {params.group}_temp2_report.tsv {params.group}_number_bins.tsv > {params.group}_temp3_report.tsv

        ls -l 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam | wc -l > {params.group}_n_samples.tsv

        nsamples=$( cat {params.group}_n_samples.tsv )
        nsamples1=$(( nsamples + 1 ))
        echo $nsamples1
        for sample in `seq 2 $nsamples1`;
            do cut -f"$sample" 3_Outputs/6_CoverM/{params.group}_assembly_coverM.txt | sed -n 3p >> {params.group}_relabun.tsv;
        done

        paste {params.group}_temp3_report.tsv {params.group}_relabun.tsv > {params.group}_temp4_report.tsv

        #Combine them into the final assembly report
        cat headers.tsv {params.group}_temp4_report.tsv > {output}

        #Clean up
#        rm *.tsv
        """