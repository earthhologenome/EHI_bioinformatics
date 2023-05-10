################################################################################
### Index coassembly
rule coassembly_index:
    input:
        os.path.join(config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"),
    output:
        os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta.rev.2.bt2l"
        ),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=96,
        time="02:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/index_assembly_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/index_assembly_log_{EHA}.log")
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