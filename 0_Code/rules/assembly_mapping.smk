################################################################################
### Map reads to assemblies
rule assembly_mapping:
    input:
        contigs=os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"
            ),
        index=os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta.rev.2.bt2l"
            ),
        r1=os.path.join(
            config["workdir"], "{PRB}/", "{EHI}_M_1.fq.gz"
            ),
        r2=os.path.join(
            config["workdir"], "{PRB}/", "{EHI}_M_2.fq.gz"
            )
    output:
        os.path.join(config["workdir"], "bams/" "{EHI}_{EHA}.bam")
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=48,
        time="05:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/assembly_mapping_benchmark_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/assembly_mapping_log_{EHI}_{EHA}.log")
    message:
        "Mapping {wildcards.EHI} to {wildcards.EHA} assembly using Bowtie2"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.contigs} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """