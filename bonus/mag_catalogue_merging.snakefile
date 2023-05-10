################################################################################
################################################################################
################################################################################
# Snakefile for annotating and merging MAG catalogues
# Raphael Eisenhofer 03/2023
#
################################################################################
################################################################################
################################################################################

### Setup sample wildcard:

MAGS = [os.path.relpath(fn, "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins").replace(".fa.gz", "")
            for fn in glob(f"3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz")]

print("Detected this many MAGs:")
len(MAGS)

################################################################################
### Setup the desired outputs
rule all:
    input:
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
        "conda_envs/3_dRep.yaml"
    threads:
        24
    resources:
        mem_gb=256,
        time='10:00:00'
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
        "conda_envs/3_GTDB-tk.yaml"
    threads:
        24
    resources:
        mem_gb=256,
        time='10:00:00'
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
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        bin = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/{MAG}.fa.gz"
    output:
        annotations = "3_Outputs/12_DRAM/{MAG}_annotations.tsv.gz",
        genes = "3_Outputs/12_DRAM/{MAG}_genes.fna.gz",
        genesfaa = "3_Outputs/12_DRAM/{MAG}_genes.faa.gz",
        genesgff = "3_Outputs/12_DRAM/{MAG}_genes.gff.gz",
        scaffolds = "3_Outputs/12_DRAM/{MAG}_scaffolds.fna.gz",
        gbk = "3_Outputs/12_DRAM/{MAG}.gbk.gz",
        distillate = directory("3_Outputs/12_DRAM/{MAG}_distillate")
    params:
        outdir = "3_Outputs/12_DRAM/{MAG}_annotate",
        mainout = "3_Outputs/12_DRAM",
        trnas = "3_Outputs/12_DRAM/{MAG}_trnas.tsv.gz",
        rrnas = "3_Outputs/12_DRAM/{MAG}_rrnas.tsv.gz",
    threads:
        2
    resources:
        mem_gb=24,
        time='04:00:00'
    benchmark:
        "3_Outputs/0_Logs/{MAG}_DRAM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{MAG}_DRAM.log"
    message:
        "Using DRAM to functionally annotate {wildcards.MAG}"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules
        
#        DRAM-setup.py import_config --config_loc /projects/mjolnir1/people/ncl550/0_software/20210705.dram.config

        DRAM.py annotate \
            -i {input.bin} \
            -o {params.outdir} \
            --threads {threads} \
#            --use_uniref \
            --min_contig_size 1500 

        #If statements for rrnas/trnas -- sometimes these won't be created
        if test -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}
        mv {params.outdir}/trnas.tsv {params.trnas}}
        else
        echo "trnas AND rrnas are both not present"
        fi

        if test -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/trnas.tsv {params.trnas}}
        else
        echo "only trnas found"
        fi

        if test ! -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}
        else
        echo "only rrnas found"
        fi

        if test ! -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            -o {output.distillate}
        else
        echo "neither trnas nor rrnas found"
        fi        

        pigz -p {threads} {params.outdir}/*.tsv
        pigz -p {threads} {params.outdir}/*.fna
        pigz -p {threads} {params.outdir}/*.faa
        pigz -p {threads} {params.outdir}/*.gff
        pigz -p {threads} {params.outdir}/genbank/*
        pigz -p {threads} {output.distillate}/*

        mv {params.outdir}/annotations.tsv.gz {output.annotations}
        mv {params.outdir}/scaffolds.fna.gz {output.scaffolds}
        mv {params.outdir}/genes.fna.gz {output.genes}
        mv {params.outdir}/*.faa.gz {output.genesfaa}
        mv {params.outdir}/*.gff.gz {output.genesgff}
        mv {params.outdir}/genbank/* {output.gbk}

        """
###############################################################################
## compress/store
rule compress:
    input:
        bin = expand("3_Outputs/12_DRAM/{MAG}_annotations.tsv.gz", MAG=MAG)
    output:
        tar = "3_Outputs/12_DRAM/DRAM.tar.gz",
    params:
        mainout = "3_Outputs/12_DRAM/"
    conda:
        "conda_envs/1_Preprocess_QC.yaml"
    threads:
        2
    resources:
        mem_gb=32,
        time='04:00:00'
    message:
        "Tarballing outputs"
    shell:
        """
        tar -cvf {output.tar} {params.mainout}
        rmdir {params.mainout}*_annotate/genbank
        rmdir {params.mainout}*_annotate
        
        """