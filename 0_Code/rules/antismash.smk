################################################################################
### Run ANTISMASH on refined bins
################################################################################
rule antismash:
    input:
        os.path.join(
            config["magdir"],
            "{MAG}.fa.gz"
        )
    output:
        annotations = os.path.join(config["magdir"], "{MAG}_anno.tsv.gz"),
        distillate = temp(directory(os.path.join(config["magdir"], "{MAG}_distillate"))),
        product = os.path.join(config["magdir"], "{MAG}_dist.tsv.gz"),
        gbk = os.path.join(config["magdir"], "{MAG}.gbk.gz")
    params:
        outdir=os.path.join(config["magdir"], "{MAG}_annotate"),
        trnas=os.path.join(config["magdir"], "{MAG}_trnas.tsv"),
        rrnas=os.path.join(config["magdir"], "{MAG}_rrnas.tsv"), 
    # conda:
    #     f"{config['codedir']}/conda_envs/DRAM.yaml"
    threads:
        2
    resources:
        mem_gb=24,
        time='03:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/DRAM_benchmark_{MAG}.tsv")
    log:
        os.path.join(config["logdir"] + "/DRAM_log_{MAG}.log")
    message:
        "Using DRAM to functionally annotate {wildcards.MAG}"
    shell:
        """
        #Loading DRAM from our custom DRAM build:
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules

            DRAM.py annotate \
                -i {input} \
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

            #compress, clean
            pigz -p {threads} {params.outdir}/annotations.tsv
            pigz -p {threads} {output.distillate}/product.tsv
            pigz -p {threads} {params.outdir}/genbank/*.gbk

            mv {params.outdir}/annotations.tsv.gz {output.annotations}
            mv {output.distillate}/product.tsv.gz {output.product}
            mv {params.outdir}/genbank/*.gbk.gz {output.gbk}
            rm -r {params.outdir}
        """