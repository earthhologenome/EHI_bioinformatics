
### Setup sample wildcard:
import os
from glob import glob

MAG = [os.path.basename(fn).replace(".fa.gz", "")
            for fn in glob(f"3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz")]

print("Detected the following MAGs:")
print(MAG)

rule all:
    input:
        expand("3_Outputs/12_DRAM/{MAG}_annotations.tsv.gz", MAG=MAG)

################################################################################
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        bin = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/{MAG}.fa.gz"
    output:
        annotation = "3_Outputs/12_DRAM/{MAG}_annotations.tsv",
    params:
        outdir = "3_Outputs/12_DRAM/{MAG}_annotate",
        mainout = "3_Outputs/12_DRAM"
    threads:
        2
    resources:
        mem_gb=32,
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

        mv {params.outdir}/* {params.mainout}

        """
################################################################################
### compress/store
rule compress:
    input:
        bin = expand("3_Outputs/12_DRAM/{MAG}_annotations.tsv")
    output:
        annotation = "3_Outputs/12_DRAM/{MAG}_annotations.tsv.gz",
    params:
        mainout = "3_Outputs/12_DRAM"
    conda:
        "conda_envs/1_Preprocess_QC.yaml"
    threads:
        32
    resources:
        mem_gb=128,
        time='04:00:00'
    message:
        "Compressing outputs"
    shell:
        """
        pigz -p {threads} {params.mainout}/*.tsv
        pigz -p {threads} {params.mainout}/*.fna
        pigz -p {threads} {params.mainout}/*.faa
        pigz -p {threads} {params.mainout}/*.gff
        pigz -p {threads} {params.mainout}/genbank/*

        tar -cvf {params.mainout}/DRAM.tar.gz {params.mainout}

        """