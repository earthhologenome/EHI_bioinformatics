################################################################################
### Bin contigs using metaWRAP's binning module
rule vamb_multisplit:
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
            config["workdir"], "{EHA}_binning/vamb_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "{EHA}_vamb")
    threads: 16
    # conda:
    #     f"{config['codedir']}/conda_envs/vamb.yaml"
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/vamb_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/vamb_log_{EHA}.log")
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

            # rename bins from .fna to .fa
            for i in {params.outdir}/bins/*.fna;
                do mv $i ${{i/.fna/.fa}};
            done

            # Create output for the next rule
            touch {output}
        """