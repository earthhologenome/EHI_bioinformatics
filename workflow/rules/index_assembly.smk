################################################################################
### Index assemblies
rule assembly_index:
    input:
        contigs=os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"),
    output:
        os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta.rev.2.bt2l"
        ),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=32,
        time="00:10:00",
    benchmark:
        os.path.join(config["logdir"] + "/index_assembly_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/index_assembly_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Indexing {wildcards.EHA} assembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input.contigs} {input.contigs} \
        &> {log}
        """