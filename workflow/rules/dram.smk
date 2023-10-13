################################################################################
### Run DRAM on refined bins
################################################################################
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        mag=os.path.join(
            config["magdir"], "{MAG}.gz"
        ),
        downloaded=os.path.join(
            config["magdir"],
            "mags_downloaded"
        )
    output:
        annotations = os.path.join(config["magdir"], "{MAG}_anno.tsv.gz"),
        product = os.path.join(config["magdir"], "{MAG}_kegg.tsv.gz"),
        gbk = os.path.join(config["magdir"], "{MAG}.gbk.gz"),
    params:
        outdir=os.path.join(config["magdir"], "{MAG}_annotate"),
        trnas=os.path.join(config["magdir"], "{MAG}_trnas.tsv"),
        rrnas=os.path.join(config["magdir"], "{MAG}_rrnas.tsv"),
        distillate=os.path.join(config["magdir"], "{MAG}_distillate")
    # conda:
    #     f"{config['codedir']}/conda_envs/DRAM.yaml"
    threads:
        2
    resources:
        mem_gb=24,
        time='04:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/DRAM_benchmark_{MAG}.tsv")
    log:
        os.path.join(config["logdir"] + "/DRAM_log_{MAG}.log")
    message:
        "Functionally annotating MAGs"
    shell:
        """
        # Remove tmp files in case job needs to be rerun
        rm -rf {params.outdir}
        rm -rf {params.distillate}

        #Loading DRAM from our custom DRAM build:
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules

            DRAM.py annotate \
                -i {input.mag} \
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
                -o {params.distillate}
            fi

            if test -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
            then
            DRAM.py distill \
                -i {params.outdir}/annotations.tsv \
                --trna_path {params.outdir}/trnas.tsv \
                -o {params.distillate}
            else
                echo "B not needed"
            fi

            if test ! -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
            then
            DRAM.py distill \
                -i {params.outdir}/annotations.tsv \
                --rrna_path {params.outdir}/rrnas.tsv \
                -o {params.distillate}
            else
                echo "C not needed"
            fi

            if test ! -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
            then
            DRAM.py distill \
                -i {params.outdir}/annotations.tsv \
                -o {params.distillate}
            else
                echo "D not needed"
            fi        

            ##get stats from annotations
            #number of gene calls
            sed '1d;' {params.outdir}/annotations.tsv | wc -l > {params.outdir}/num_genes.tsv
            sleep 5
            #number of genes with KEGG hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f9 | grep -v '^$' | wc -l > {params.outdir}/kegg_hits.tsv
            sleep 5
            #number of pfam hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f18 | grep -v '^$' | wc -l > {params.outdir}/pfam_hits.tsv
            sleep 5
            #number of cazy hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f19 | grep -v '^$' | wc -l > {params.outdir}/cazy_hits.tsv
            sleep 5
            #number of genes without annotations
            sleep 5
            bash {config[codedir]}/scripts/count_unannotated.sh {params.outdir}/annotations.tsv {params.outdir}/unannotated.tsv 2> {params.outdir}/error.log
            sleep 5

            #fetch only KEGG modules from product
            cut -f1-381 {params.distillate}/product.tsv > {params.distillate}/kegg_product.tsv

            #compress, clean, move
            pigz -p {threads} {params.outdir}/annotations.tsv
            pigz -p {threads} {params.distillate}/kegg_product.tsv
            pigz -p {threads} {params.outdir}/genbank/*.gbk

            mv {params.outdir}/annotations.tsv.gz {output.annotations}
            mv {params.distillate}/kegg_product.tsv.gz {output.product}
            mv {params.outdir}/genbank/*.gbk.gz {output.gbk}
          
        """