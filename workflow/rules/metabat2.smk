################################################################################
### Bin contigs using metabat2
rule metabat2:
    input:
        bam=os.path.join(
            config["workdir"], "bams/", "{PRB}_{EHI}_{EHA}.bam"
            ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{PRB}_{EHI}_{EHA}_binning/metabat2_binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_binning"),
        contigsize=config["contigsize"]
    conda:
        f"{config['codedir']}/conda_envs/metabat2.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/metabat2_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/metabat2_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with metabat2"
    shell:
        """
            # summarise contig depths
            jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/metabat_depth.txt {input.bam}

            {params.outdir}/metabat2_bins/

            # Run metabat2
            metabat2 \
            -i {input.contigs} \
            -a {params.outdir}/metabat_depth.txt \
            -o {params.outdir/metabat2_bins/bin \
            -m 1500 \
            -t {threads} \
            --unbinned

            # Create output for the next rule
            touch {output}

        """