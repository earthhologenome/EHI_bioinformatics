################################################################################
### Run DRAM on refined bins
################################################################################
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        os.path.join(
            config["magdir"],
            "{MAG}.fa.gz"
        )
    output:
        annotations = os.path.join(config["magdir"], "{MAG}_anno.tsv.gz"),
        genes = temp(os.path.join(config["magdir"], "{MAG}_genes.fna.gz")),
        genesfaa = temp(os.path.join(config["magdir"], "{MAG}_genes.faa.gz")),
        genesgff = temp(os.path.join(config["magdir"], "{MAG}_genes.gff.gz")),
        scaffolds = temp(os.path.join(config["magdir"], "{MAG}_scaffolds.fna.gz")),
        gbk = temp(os.path.join(config["magdir"], "{MAG}.gbk.gz")),
        distillate = directory(os.path.join(config["magdir"], "{MAG}_distillate")),
        product = os.path.join(config["magdir"], "{MAG}_dist.tsv.gz")
    params:
        outdir=os.path.join(config["magdir"], "{MAG}_annotate"),
        trnas=os.path.join(config["magdir"], "{MAG}_trnas.tsv.gz"),
        rrnas=os.path.join(config["magdir"], "{MAG}_rrnas.tsv.gz"), 
    conda:
        f"{config['codedir']}/conda_envs/DRAM.yaml"
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
        "Using DRAM to functionally annotate {wildcards.MAG}"
    shell:
        """
            DRAM.py annotate \
                -i {input.mags} \
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

            #calculate number of contigs and genes for later
            grep '>' {params.outdir}/scaffolds.fna | wc -l > {config[workdir}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}/DRAM/{input.mags}_ncontigs.txt
            grep '>' {params.outdir}/genes.fna | wc -l > {config[workdir}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}/DRAM/{input.mags}_ngenes.txt

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
            mv {output.distillate}/product.tsv.gz {output.product}

            rm {output.distillate}/*

            touch {output.complete}
        """