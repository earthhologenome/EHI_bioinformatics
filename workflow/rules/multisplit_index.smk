################################################################################
### Index multisplit contigs
rule multisplit_index:
    input:
        os.path.join(config["workdir"], "{EHA}_combined_contigs.fasta.gz")
    output:
        os.path.join(config["workdir"], "{EHA}_combined_contigs.fasta.gz.rev.2.bt2l")
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time="02:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/{EHA}_index_multisplit_benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/{EHA}_index_multisplit_log.log")
    message:
        "Indexing {wildcards.EHA} assembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {input} \
        &> {log}
        """