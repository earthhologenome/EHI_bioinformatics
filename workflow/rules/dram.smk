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
        distillate = directory(os.path.join(config["magdir"], "{MAG}_distillate"))
    params:
        outdir=os.path.join(config["magdir"], "{MAG}_annotate"),
        trnas=os.path.join(config["magdir"], "{MAG}_trnas.tsv"),
        rrnas=os.path.join(config["magdir"], "{MAG}_rrnas.tsv")
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
        rm -rf {output.distillate}

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
                -o {output.distillate}
            else
            echo "trnas AND rrnas are both not present"
            fi

            if test -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
            then
            DRAM.py distill \
                -i {params.outdir}/annotations.tsv \
                --trna_path {params.outdir}/trnas.tsv \
                -o {output.distillate}
            else
            echo "only trnas found"
            fi

            if test ! -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
            then
            DRAM.py distill \
                -i {params.outdir}/annotations.tsv \
                --rrna_path {params.outdir}/rrnas.tsv \
                -o {output.distillate}
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

            ##get stats from annotations
            #number of gene calls
            sed '1d;' {params.outdir}/annotations.tsv | wc -l > {params.outdir}/num_genes.tsv
            #number of genes with KEGG hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f9 | grep -v '^$' | wc -l > {params.outdir}/kegg_hits.tsv
            #number of pfam hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f18 | grep -v '^$' | wc -l > {params.outdir}/pfam_hits.tsv
            #number of cazy hits
            sed '1d;' {params.outdir}/annotations.tsv | cut -f19 | grep -v '^$' | wc -l > {params.outdir}/cazy_hits.tsv
            #number of genes without annotations
            awk -F'\t' 'BEGIN {{ count = 0; }} {{ empty = 1; for (i = 9; i <= 19; i++) {{ if ($i != "") {{ empty = 0; break; }} }} if (empty == 1) {{ count++; }} }} END {{ print count; }}' {params.outdir}/annotations.tsv > {params.outdir}/unannotated.tsv



            #fetch only KEGG modules from product (keep DRAM distillate)
            cut -f1-381 {output.distillate}/product.tsv > {output.distillate}/kegg_product.tsv
#            cut -f1,382-466 {output.distillate}/product.tsv > {output.distillate}/dram_product.tsv

            #compress, clean, move
            pigz -p {threads} {params.outdir}/annotations.tsv
            pigz -p {threads} {output.distillate}/kegg_product.tsv
            pigz -p {threads} {params.outdir}/genbank/*.gbk

            mv {params.outdir}/annotations.tsv.gz {output.annotations}
            mv {output.distillate}/kegg_product.tsv.gz {output.product}
            mv {params.outdir}/genbank/*.gbk.gz {output.gbk}
          
        """