################################################################################
################################################################################
################################################################################
# Snakefile for MAG dereplication, taxonomic and functional annotation.
# Raphael Eisenhofer 02/2022
#
################################################################################
################################################################################
################################################################################

configfile: "0_Code/3_dRep_config.yaml"

### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_M_1.fastq.gz", "")
            for fn in glob(f"2_Reads/4_Host_removed/*_M_1.fastq.gz")]

MAGS = [os.path.relpath(fn, "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins").replace(".fa.gz", "")
            for fn in glob(f"3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz")]

print("Detected the following samples:")
print(SAMPLE)

print("Detected this many MAGs:")
len(MAGS)
################################################################################
### Setup the desired outputs
rule all:
    input:
#        expand("3_Outputs/11_DRAM/Distillate/DRAM_product.html")
        "3_Outputs/10_Final_tables/unfiltered_count_table.txt"
################################################################################
### Dereplicate refined bins using dRep
rule dereplication:
    input:
        bins = "3_Outputs/5_Refined_Bins/"
    output:
        drep = "3_Outputs/7_Dereplication/figures/Primary_clustering_dendrogram.pdf"
    params:
        ANI = expand("{ANI}", ANI=config['ANI']),
        workdir = "3_Outputs/7_Dereplication/"
    conda:
        "3_dRep.yaml"
    threads:
        40
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/dRep.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/dRep.log"
    message:
        "Dereplicating bins that are > {params.ANI} percent indentical"
    shell:
        """
        # Parse/collate metawrap stats files for compatibility with dRep genomeinfo:
        echo -e "genome,completeness,contamination" > {input.bins}/header.txt
        for i in {input.bins}/*.stats;
            do sed '1d;' $i | cut -f 1,2,3 --output-delimiter ',' >> {input.bins}/bin_info.txt;
                done
        sed -i'' 's@^@{input.bins}/All_metawrap_70_10_bins/@g' {input.bins}/bin_info.txt
        sed -i'' 's/,/.fa.gz,/' {input.bins}/bin_info.txt
        cat {input.bins}/header.txt {input.bins}/bin_info.txt > {input.bins}/genome_info.csv
        rm {input.bins}/*.txt

# Decompress bins (dRep can't handle .gz input) -- learned this the hard way!
# gunzip {input.bins}/All_metawrap_70_10_bins/*.fa.gz

        # Run dRep
            dRep dereplicate \
                {params.workdir} \
                -p {threads} \
                -comp 70 \
                -sa {params.ANI} \
                -g {input.bins}/All_metawrap_70_10_bins/*.fa.gz \
                --genomeInfo {input.bins}/genome_info.csv
                2> {log}

# Rename output, compress bins
# pigz -p {threads} {input.bins}/All_metawrap_70_10_bins/*.fa
# pigz -p {threads} {params.workdir}/dereplicated_genomes/*.fa
        """
################################################################################
### Annotate dereplicated MAGs with gtdb-tk taxonomy:
rule gtdbtk:
    input:
        "3_Outputs/7_Dereplication/figures/Primary_clustering_dendrogram.pdf"
    output:
        "3_Outputs/8_GTDB-tk/classify/gtdbtk.bac120.summary.tsv"
    params:
        GTDB_data = expand("{GTDB_data}", GTDB_data=config['GTDB_data']),
        outdir = "3_Outputs/8_GTDB-tk/",
        bins = "3_Outputs/7_Dereplication/dereplicated_genomes"
    conda:
        "3_GTDB-tk.yaml"
    threads:
        40
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/gtdbtk.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/gtdbtk.log"
    message:
        "Running gtdb-tk to taxonomically annotate MAGs"
    shell:
        """
        # Specify path to reference data:
        export GTDBTK_DATA_PATH={params.GTDB_data}

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.bins} \
        --extension "gz" \
        --out_dir {params.outdir} \
        --cpus {threads}

        # Create a merged summary output for DRAM:
        if [ -s "{params.outdir}/classify/gtdbtk.ar122.summary.tsv" ]
        then
        sed '1d;' {params.outdir}/classify/gtdbtk.ar122.summary.tsv > {params.outdir}/ar122.tsv
        cat {output} {params.outdir}/ar122.tsv > {params.outdir}/gtdbtk_combined_summary.tsv
        rm {params.outdir}/ar122.tsv

        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {output} > {params.outdir}/gtdbtk_combined_summary.tsv
        fi
        """
################################################################################
### Index the MAG catalogue
rule Coassembly_index:
    input:
        "3_Outputs/8_GTDB-tk/classify/gtdbtk.bac120.summary.tsv"
    output:
        "3_Outputs/9_MAG_catalogue_mapping/MAGs.fa.gz.rev.2.bt2l"
    params:
        MAGs = "3_Outputs/7_Dereplication/dereplicated_genomes",
        catted_MAGs = "3_Outputs/9_MAG_catalogue_mapping/MAGs.fa.gz",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/MAG_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/MAG_indexing.log"
    message:
        "Indexing MAG catalogue using Bowtie2"
    shell:
        """
        # Rename dereplicated MAG headers for CoverM compatibility
        for i in {params.MAGs}/*.fa.gz;
            do rename.sh \
                in=$i \
                out=${{i/.fa.gz/_renamed.fa.gz}} \
                zl=9 \
                prefix=$(basename ${{i/.fa.gz/^}});
            done

        # Concatenate the dereplicated MAGs into a single file
        cat {params.MAGs}/*_renamed.fa.gz > {params.catted_MAGs}

        # Index the MAG catalogue
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.catted_MAGs} {params.catted_MAGs} \
        &> {log}

        # Clean up
        rm {params.MAGs}/*_renamed.fa.gz
        """
################################################################################
### Map the preprocessed reads to the dereplicated MAG catalogue
rule MAG_catalogue_mapping:
    input:
        index = "3_Outputs/9_MAG_catalogue_mapping/MAGs.fa.gz.rev.2.bt2l",
        r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq.gz",
        r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq.gz"
    output:
        "3_Outputs/9_MAG_catalogue_mapping/BAMs/{sample}.bam"
    params:
        BAMs = "3_Outputs/9_MAG_catalogue_mapping/BAMs",
        MAGs = "3_Outputs/9_MAG_catalogue_mapping/MAGs.fa.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        20
    resources:
        mem_gb=90
    benchmark:
        "3_Outputs/0_Logs/{sample}_MAG_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_MAG_mapping.log"
    message:
        "Mapping {wildcards.sample} to the dereplicated MAG catalogue using Bowtie2"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.MAGs} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """
################################################################################
### Create the final count table using CoverM
rule coverM_assembly:
    input:
        expand("3_Outputs/9_MAG_catalogue_mapping/BAMs/{sample}.bam", sample=SAMPLE)
    output:
        "3_Outputs/10_Final_tables/unfiltered_count_table.txt"
    params:
        BAMs = "3_Outputs/9_MAG_catalogue_mapping/BAMs",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    resources:
        mem_gb=180
    benchmark:
        "3_Outputs/0_Logs/coverm.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/coverm.log"
    message:
        "Creating the count table with CoverM"
    shell:
        """
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_fraction length \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
################################################################################
### Create genome-scale metabolic models (GEMs)
### Adapted from metaGEM: 
### https://github.com/franciscozorrilla/metaGEM/wiki/Reconstruct-&-evaluate-genome-scale-metabolic-models-with-CarveMe-and-memote
# rule extractProteinBins:
#     message:
#         "Extract ORF annotated protein fasta files for each bin from reassembly checkm files."
#     shell:
#         """
#         mkdir -p 3_Outputs/12_protein_bins

#         echo -e "Begin moving and renaming ORF annotated protein fasta bins from dereplicated_genomes/ to protein_bins/ ... \n"
#         for bin in 3_Outputs/7_Dereplication/dereplicated_genomes/*.fa.gz;
#             do var=$(echo $bin/genes.faa | sed 's|reassembled_bins/||g'|sed 's|/reassembled_bins.checkm/bins||'|sed 's|/genes||g'|sed 's|/|_|g'|sed 's/permissive/p/g'|sed 's/orig/o/g'|sed 's/strict/s/g');
#             cp $bin/*.faa {config[path][root]}/{config[folder][proteinBins]}/$var;
#         done
#         """


################################################################################
### Setup DRAM groups to split into multiple jobs for efficiency.
### We can use snakemake checkpoints for this :)
#     input:
#         "samples/{sample}.txt"
#     output:
#         clusters=directory("clustering/{sample}")
#     shell:
#         """
#         "mkdir clustering/{wildcards.sample}; "
#         "for i in 1 2 3; do echo $i > clustering/{wildcards.sample}/$i.txt; done"
#
#
#         ## Split bins into groups to parallelize DRAM:
#         # How many bins?
#         count=$(find {params.bins}/ -name '*.fa.gz' -type f|wc -l)
#         # How many bins per group (using 5 groups)?
#         groupsize=$(((count +4) / 5))
#         # Move bins into separate group folders:
#         for group in `seq 1 5`;
#             do mkdir -p {params.workdir}/"group$group";
#             find {params.bins} -type f | head -n $groupsize |
#             xargs -i mv "{{}}" {params.workdir}/"group$group"; done
#         """
# ################################################################################
# ### input function for rule aggregate, return paths to all files produced by the
# ### checkpoint 'somestep'
# # an intermediate rule
# rule intermediate:
#     input:
#         "clustering/{sample}/{i}.txt"
#     output:
#         "post/{sample}/{i}.txt"
#     shell:
#         "cp {input} {output}"
#
#
# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
#     return expand("post/{sample}/{i}.txt",
#            sample=wildcards.sample,
#            i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)
#
#
# # an aggregation over all produced clusters
# rule aggregate:
#     input:
#         aggregate_input
#     output:
#         "aggregated/{sample}.txt"
#     shell:
#         "cat {input} > {output}"
# ################################################################################
# ### Run parallel DRAM annotate jobs.
# rule DRAM_annotate:
#     input:
#         aggregate_input
#     output:
#         "3_Outputs/11_DRAM/Distillate/product.html"
#     params:
#         DRAM_data = expand("{DRAM_data}", DRAM_data=config['DRAM_data']),
#         bins = "3_Outputs/7_Dereplication/dereplicated_genomes",
#         workdir = "3_Outputs/11_DRAM/",
#         gtdbtax = "3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv",
#         checkm_stats =
#     conda:
#         "3_DRAM.yaml"
#     threads:
#         40
#     resources:
#         mem_gb=180
#     benchmark:
#         "3_Outputs/0_Logs/DRAM.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/DRAM.log"
#     message:
#         "Functionally annotating MAGs using DRAM"
#     shell:
#         """
#         # Create checkm tsv for input to DRAM:
#         echo -e "Bin Id\tCompleteness\tContamination" > {params.workdir}/header.txt
#         sed '1d;' 3_Outputs/5_Refined_Bins/dRep_groups/{wildcards.group}/genome_info.csv |
#         tr ',' '\t' > {params.workdir}/bininfo.txt
#         cat {params.workdir}/header.txt {params.workdir}/bininfo.txt > {params.workdir}checkm.tsv
#
#         # Clean up
#         rm {params.workdir}/*.txt
#
#         # Run DRAM-annotate:
#         DRAM.py annotate \
#         -i '$group/*.fa.gz' \
#         --gtdb_taxonomy {params.gtdbtax} \
#         --checkm_quality {params.workdir}checkm.tsv \
#         --threads 8 \
#         --min_contig_size 2500 \
#         -o "$group"_DRAM
#         """
# ################################################################################
# ###
#
#
#         # Merge DRAM groups
#         DRAM.py merge_annotation \
#         -i {params.workdir}/'group*_DRAM' \
#         -o {params.workdir}/merged_DRAM
#
#         # Distill annotations:
#         DRAM.py distill \
#         -i {params.workdir}/merged_DRAM/annotations.tsv \
#         --rrna_path {params.workdir}/merged_DRAM/rrnas.tsv \
#         --trna_path {params.workdir}/merged_DRAM/trnas.tsv \
#         -o {params.workdir}/Distillate
#
#         # Rename, clean, compress:
#         for i in {params.workdir}/Distillate/*;
#             do mv $i {params.workdir}/Distillate/$(basename {wildcards.group}_"$i");
#                 done
#         pigz -p {threads} {params.workdir}/Distillate/*
#         rm -r {params.workdir}/group*_DRAM
