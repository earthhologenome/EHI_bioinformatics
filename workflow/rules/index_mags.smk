################################################################################
### Index mags
rule mag_index:
    input:
        drep=os.path.join(
            config["workdir"],
            "drep/",
            "figures/",
            config["dmb"] + "_Primary_clustering_dendrogram.pdf"
        )
    output:
        mags=os.path.join(
            config["workdir"],
            "mag_catalogue/",
            config["dmb"] + "_mags.fasta.gz"
        ),
        index=os.path.join(
            config["workdir"],
            "mag_catalogue/",
            config["dmb"] + "_mags.fasta.gz.rev.2.bt2l"
        )
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=96,
        time="02:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/index_mags_benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/index_mags_log.log")
    message:
        "Indexing mags using Bowtie2"
    shell:
        """
        # Concatenate MAGs
        cat {config[workdir]}/drep/dereplicated_genomes/*.fa.gz > {config[workdir]}/mag_catalogue/{config[dmb]}_mags.fasta.gz

        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.mags} {output.mags} \
        &> {log}
        """