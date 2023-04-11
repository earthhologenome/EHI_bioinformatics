################################################################################
### Index assemblies
rule assembly_index:
    input:
        os.path.join(config["workdir"], "{PRB}" "{EHI}" "{EHA}_contigs.fasta"),
    output:
        os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_contigs.fasta.rev.2.bt2l"
        ),
    conda:
        f"{config['codedir']}/conda_envs/2_Assembly_Binning.yaml"
    threads: 16
    resources:
        mem_gb=96,
        time="02:00:00",
    benchmark:
        "{{config['logdir']}}/index_assembly_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/index_assembly_log_{PRB}_{EHI}_{EHA}.log",
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