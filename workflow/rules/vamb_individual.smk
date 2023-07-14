################################################################################
### Bin contigs using metaWRAP's binning module
rule vamb_individual:
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
        outdir=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_vamb")
    threads: 16
    # conda:
    #     f"{config['codedir']}/conda_envs/vamb.yaml"
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/vamb_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/vamb_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with vamb"
    shell:
        """
            conda activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/vamb4.1.3
            # Run vamb
            vamb \
            --outdir {params.outdir} \
            --fasta {input.contigs} \
            --bamfiles {input.bam} \
            --minfasta 200000

            # Create output for the next rule
            touch {output}
        """