################################################################################
################################################################################
################################################################################
# Snakefile for MAG dereplication, taxonomic and functional annotation.
# Raphael Eisenhofer 02/2022
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:
import os
from glob import glob

GROUP = [os.path.basename(dir)
         for dir in glob(f"3_Outputs/5_Refined_Bins/dRep_groups/*")]

MAGS = [os.path.relpath(fn, "3_Outputs/5_Refined_Bins/").replace(".fa.gz", "")
            for group in GROUP
            for fn in glob(f"3_Outputs/5_Refined_Bins/{group}/bins/*.fa.gz")]

print("Detected these sample groups:")
print(GROUP)
print("Detected this many MAGs:")
len(MAGS)
################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("3_Outputs/11_DRAM/{group}/Distillate/{group}_product.html", group=GROUP)
################################################################################
### Dereplicate refined bins using dRep
rule dereplication:
    input:
        bins = "3_Outputs/5_Refined_Bins/dRep_groups/{group}"
    output:
        drep = "3_Outputs/7_Dereplication/{group}/figures/{group}_Primary_clustering_dendrogram.pdf"
    params:
        ANI = expand("{ANI}", ANI=config['ANI']),
        workdir = "3_Outputs/7_Dereplication/{group}"
    conda:
        "3_dRep.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_dRep.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_dRep.log"
    message:
        "Dereplicating bins for {wildcards.group} that are > {params.ANI} percent indentical"
    shell:
        """
        # Parse/collate metawrap stats files for compatibility with dRep genomeinfo:
        echo -e "genome,completeness,contamination" > {input.bins}/header.txt
        for i in {input.bins}/*.stats;
            do sed '1d;' $i | cut -f 1,2,3 --output-delimiter ',' >> {input.bins}/bin_info.txt;
                done
        sed -i'' 's@^@{input.bins}/bins/@g' {input.bins}/bin_info.txt
        sed -i'' 's/,/.fa,/' {input.bins}/bin_info.txt
        cat {input.bins}/header.txt {input.bins}/bin_info.txt > {input.bins}/genome_info.csv
        rm {input.bins}/*.txt

        # Decompress bins (dRep can't handle .gz input) -- learned this the hard way!
        gunzip {input.bins}/bins/*.fa.gz

        # Run dRep
            dRep dereplicate \
                {params.workdir} \
                -p {threads} \
                -comp 70 \
                -sa {params.ANI} \
                -g {input.bins}/bins/*.fa \
                --genomeInfo {input.bins}/genome_info.csv
                2> {log}

        # Rename output, compress bins
        for i in {params.workdir}/figures/*; do
            mv $i {params.workdir}/figures/{wildcards.group}_$(basename $i);
                done
        pigz -p {threads} {input.bins}/bins/*.fa
        pigz -p {threads} {params.workdir}/dereplicated_genomes/*.fa
        """
################################################################################
### Annotate dereplicated MAGs with gtdb-tk taxonomy:
rule gtdbtk:
    input:
        "3_Outputs/7_Dereplication/{group}/figures/{group}_Primary_clustering_dendrogram.pdf"
    output:
         "3_Outputs/8_GTDB-tk/{group}/classify/{group}.bac120.summary.tsv"
    params:
        outdir = "3_Outputs/8_GTDB-tk/{group}",
        bins = "3_Outputs/7_Dereplication/{group}/dereplicated_genomes"
    conda:
        "3_GTDB-tk.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_gtdbtk.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_gtdbtk.log"
    message:
        "Running gtdb-tk on {wildcards.group} MAGs"
    shell:
        """
        # Specify path to reference data:
        GTDBTK_DATA_PATH=/home/projects/ku-cbd/people/rapeis/0_DBs/release207_v2

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.bins} \
        --extension "gz" \
        --out_dir {params.outdir} \
        --cpus {threads} \
        --scratch_dir {params.outdir} \
        --prefix {wildcards.group}

        # Create a merged summary output for DRAM:
        if [ -s "{params.outdir}/classify/{wildcards.group}.ar122.summary.tsv" ]
        then
        sed '1d;' {params.outdir}/classify/{wildcards.group}.ar122.summary.tsv > {params.outdir}/ar122.tsv
        cat {output} {params.outdir}/ar122.tsv > {params.outdir}/gtdbtk_combined_summary.tsv
        rm {params.outdir}/ar122.tsv

        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {output} > {params.outdir}/{wildcards.group}_combined_summary.tsv
        fi
        """
################################################################################
### Index the MAG catalogue
rule Coassembly_index:
    input:
        "3_Outputs/8_GTDB-tk/{group}/classify/{group}.bac120.summary.tsv"
    output:
        "3_Outputs/9_MAG_catalogue_mapping/{group}/{group}_MAGs.fa.gz.rev.2.bt2l"
    params:
        MAGs = "3_Outputs/7_Dereplication/{group}/dereplicated_genomes",
        catted_MAGs = "3_Outputs/9_MAG_catalogue_mapping/{group}/{group}_MAGs.fa.gz",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_MAG_indexing.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_MAG_indexing.log"
    message:
        "Indexing {wildcards.group} MAG catalogue using Bowtie2"
    shell:
        """
        # Rename dereplicated MAG headers for CoverM compatibility
        for i in {params.MAGs}/*.fa.gz;
            do rename.sh \
                in=$i \
                out=${i/.fa.gz/_renamed.fa.gz} \
                zl=9 \
                prefix=$(basename ${i/.fa.gz/-});
            done

        # Concatenate the dereplicated MAGs into a single file
        cat {params.MAGs}/*_renamed.fa.gz > {params.catted_MAGs}

        # Index the MAG catalogue
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.catted_MAGs} {params.catted_MAGs} \
        &> {log}
        """
################################################################################
### Map the preprocessed reads to the dereplicated MAG catalogue
rule MAG_catalogue_mapping:
    input:
        "3_Outputs/9_MAG_catalogue_mapping/{group}/{group}_MAGs.fa.gz.rev.2.bt2l"
    output:
        "3_Outputs/9_MAG_catalogue_mapping/{group}/BAMs/Done.txt"
    params:
        reads = "2_Reads/3_Host_removed/{group}",
        BAMs = "3_Outputs/9_MAG_catalogue_mapping/{group}/BAMs",
        MAGs = "3_Outputs/9_MAG_catalogue_mapping/{group}/{group}_MAGs.fa.gz"
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_MAG_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_MAG_mapping.log"
    message:
        "Mapping {wildcards.group} samples to MAG their MAG catalogue using Bowtie2"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        for fq1 in {params.reads}/*_1.fastq.gz; do \
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.MAGs} \
            -1 $fq1 \
            -2 ${{fq1/_1.fastq.gz/_2.fastq.gz}} \
        | samtools sort -@ {threads} -o {params.BAMs}/$(basename ${{fq1/_1.fastq.gz/.bam}}); done

        #Create output file for snakemake
        echo "hi" > {output}
        """
################################################################################
### Create the final count table using CoverM
rule coverM_assembly:
    input:
        "3_Outputs/9_MAG_catalogue_mapping/{group}/BAMs/Done.txt"
    output:
        "3_Outputs/10_Final_tables/{group}_final_count_table.txt"
    params:
        BAMs = "3_Outputs/9_MAG_catalogue_mapping/{group}/BAMs",
    conda:
        "2_Assembly_Binning.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_coverm.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coverm.log"
    message:
        "Creating final count table for {wildcards.group} with CoverM"
    shell:
        """
        coverm genome \
            -b {params.BAMs}/*.bam \
            -s - \
            -m count covered_fraction length \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """
################################################################################
### Functionally annotate/distil MAGs with DRAM
rule DRAM_annotate:
    input:
        "3_Outputs/10_Final_tables/{group}_final_count_table.txt"
    output:
        "3_Outputs/11_DRAM/{group}/Distillate/{group}_product.html"
    params:
        bins = "3_Outputs/7_Dereplication/{group}/dereplicated_genomes",
        workdir = "3_Outputs/11_DRAM/{group}",
        gtdbtax = "3_Outputs/8_GTDB-tk/{group}/gtdbtk_combined_summary.tsv",
    conda:
        "3_DRAM.yaml"
    threads:
        40
    benchmark:
        "3_Outputs/0_Logs/{group}_DRAM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_DRAM.log"
    message:
        "Functionally annotating {wildcards.group} MAGs using DRAM"
    shell:
        """
        ## Split bins into groups to parallelize DRAM:
        # How many bins?
        count=$(find {params.bins}/ -name '*.fa.gz' -type f|wc -l)
        # How many bins per group (using 5 groups)?
        groupsize=$(((count +4) / 5))
        # Move bins into separate group folders:
        for group in `seq 1 5`;
            do mkdir -p {params.workdir}/"group$group";
            find {params.bins} -type f | head -n $groupsize |
            xargs -i mv "{{}}" {params.workdir}/"group$group"; done

        # Create checkm tsv for input to DRAM:
        echo -e "Bin Id\tCompleteness\tContamination" > {params.workdir}/header.txt
        sed '1d;' 3_Outputs/5_Refined_Bins/dRep_groups/{wildcards.group}/genome_info.csv |
        tr ',' '\t' > {params.workdir}/bininfo.txt
        cat {params.workdir}/header.txt {params.workdir}/bininfo.txt > {params.workdir}checkm.tsv

        # Clean up
        rm {params.workdir}/*.txt

        # Run DRAM-annotate:
        for group in {params.workdir}/group*;
        do echo DRAM.py annotate \
        -i '$group/*.fa.gz' \
        --gtdb_taxonomy {params.gtdbtax} \
        --checkm_quality {params.workdir}checkm.tsv \
        --threads 8 \
        --min_contig_size 2500 \
        -o "$group"_DRAM;
            done | parallel -j 5

        # Merge DRAM groups
        DRAM.py merge_annotation \
        -i {params.workdir}/'group*_DRAM' \
        -o {params.workdir}/merged_DRAM

        # Distill annotations:
        DRAM.py distill \
        -i {params.workdir}/merged_DRAM/annotations.tsv \
        --rrna_path {params.workdir}/merged_DRAM/rrnas.tsv \
        --trna_path {params.workdir}/merged_DRAM/trnas.tsv \
        -o {params.workdir}/Distillate

        # Rename, clean, compress:
        for i in {params.workdir}/Distillate/*;
            do mv $i {params.workdir}/Distillate/$(basename {wildcards.group}_"$i");
                done
        pigz -p {threads} {params.workdir}/Distillate/*
        rm -r {params.workdir}/group*_DRAM
        """
################################################################################
###
