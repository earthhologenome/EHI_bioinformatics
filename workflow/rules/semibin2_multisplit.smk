################################################################################
### Bin contigs using metaWRAP's binning module
rule semibin2_individual:
    input:
        bam=expand(
                os.path.join(
                config["workdir"], "bams/", "{combo[0]}_{combo[1]}_{combo[2]}.bam"
                ), combo=valid_combinations
            ),
        contigs=os.path.join(
            config["workdir"], 
            "{EHA}_combined_contigs.fasta.gz"
            )
    output:
        os.path.join(
            config["workdir"], "{EHA}_binning/semibin2_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{EHA}_semibin2")
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    # conda:
    #     f"{config['codedir']}/conda_envs/semibin2.yaml"
    benchmark:
        os.path.join(config["logdir"] + "/semibin2_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/semibin2_log_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with semibin2"
    shell:
        """
            conda activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/semibin_1.5.1
            # Run semibin2
            SemiBin2 multi_easy_bin \
                    --self-supervised \
                    -i {input.contigs} \
                    -b {input.bam} \
                    -o {params.outdir} \
                    -p {threads} -t {threads} \
                    --separator :
            
            #gunzip bins for metawrap refinement
            gunzip {params.outdir}/output_bins/*

            # Create output for the next rule
            touch {output}

        """