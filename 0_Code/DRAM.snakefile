
### Setup sample wildcard:
import os
from glob import glob

MAG = [os.path.basename(fn).replace(".fa.gz", "")
            for fn in glob(f"3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz")]

print("Detected the following MAGs:")
print(MAG)

rule all:
    input:
        expand("3_Outputs/12_DRAM/{MAG}/{MAG}_annotations.tsv", MAG=MAG)

################################################################################
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        bin = "3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/{MAG}.fa.gz"
    output:
        annotation = "3_Outputs/12_DRAM/{MAG}/{MAG}_annotations.tsv.gz",
    params:
        outdir = "3_Outputs/12_DRAM/{MAG}"
    # conda:
    #     "conda_envs/3_DRAM.yaml"
    threads:
        2
    resources:
        mem_gb=32,
        time='03:00:00'
    benchmark:
        "3_Outputs/0_Logs/{MAG}_DRAM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{MAG}_DRAM.log"
    message:
        "Using DRAM to functionally annotate bin"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules
        
#        DRAM-setup.py import_config --config_loc /projects/mjolnir1/people/ncl550/0_software/20210705.dram.config


        DRAM.py annotate \
            -i {input.bin} \
            -o {output.outdir} \
            --threads {threads} \
            --min_contig_size 1500 

        pigz -p {threads} {params.outdir}/*

        for i in {params.outdir}/*; do mv $i {params.outdir}/{MAG}_$(basename $i); done

        """