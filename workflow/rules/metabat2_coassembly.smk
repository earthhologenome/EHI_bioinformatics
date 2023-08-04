################################################################################
### Bin contigs using metabat2
rule metabat2:
    input:
        bam=expand(os.path.join(
            config["workdir"], "bams/", "{combo[0]}_{combo[1]}_{combo[2]}.bam"
            ),
            combo=valid_combinations
        ),
        contigs=os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{EHA}_binning/metabat2_binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{EHA}_binning"),
        contigsize=config["contigsize"]
    conda:
        f"{config['codedir']}/conda_envs/metabat2.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/binning_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/binning_log_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with metabat2"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt {params.contigsize} ]
        then
            touch {output}

        else

            # summarise contig depths
            jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/metabat_depth.txt {input.bam}

            {params.outdir}/metabat2_bins/

            # Run metabat2
            metabat2 \
            -i {input.contigs} \
            -a {params.outdir}/metabat_depth.txt \
            -o {params.outdir/metabat2_bins/bin \
            -m 1500 \
            -t {params.threads} \
            --unbinned

            # Create output for the next rule
            touch {output}

        fi
        """