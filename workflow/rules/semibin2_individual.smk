################################################################################
### Bin contigs using metaWRAP's binning module
rule semibin2_individual:
    input:
        bam=os.path.join(
            config["workdir"], "bams/", "{PRB}_{EHI}_{EHA}.bam"
            ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{PRB}_{EHI}_{EHA}_binning/binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_semibin2")
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    # conda:
    #     f"{config['codedir']}/conda_envs/semibin2.yaml"
    benchmark:
        os.path.join(config["logdir"] + "/semibin2_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/semibin2_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with semibin2"
    shell:
        """
            conda activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/semibin_1.5.1
            # Run semibin2
            SemiBin single_easy_bin \
                    --training-type self \
                    -i {input.contigs} \
                    -b {input.bam} \
                    -o {params.outdir}

            # Create output for the next rule
            touch {output}

        """