
### Setup sample wildcard:
import os
from glob import glob

MAG = [os.path.basename(fn).replace(".fa.gz", "")
            for fn in glob(f"3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz")]

print("Detected the following MAGs:")
print(MAG)

rule all:
    input:
        expand("3_Outputs/12_DRAM/{MAG}_annotations.tsv.gz", MAG=MAG),
        "3_Outputs/12_DRAM/DRAM.tar.gz"


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


        pigz -p {threads} {params.outdir}/*.tsv
        pigz -p {threads} {params.outdir}/*.fna
        pigz -p {threads} {params.outdir}/*.faa
        pigz -p {threads} {params.outdir}/*.gff
        pigz -p {threads} {params.outdir}/genbank/*

        #If statements for rrnas/trnas -- sometimes these won't be created
        if test -f {params.trnas}
        then 
        mv {params.outdir}/trnas.tsv.gz {output.trnas}
        else
        echo "no trnas file"
        fi

        if test -f {params.rrnas}
        then 
        mv {params.outdir}/rrnas.tsv.gz {output.rrnas}
        else
        echo "no rrnas file"
        fi

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
        32
    resources:
        mem_gb=128,
        time='04:00:00'
    message:
        "Tarballing outputs"
    shell:
        """
        tar -cvf {output.tar} {params.mainout}
        rmdir {params.mainout}/*_annotate/genbank
        rmdir {params.mainout}/*_annotate
        """