################################################################################
### Map reads to multisplit assembly
rule multisplit_mapping:
    input:
        index=os.path.join(
            config["workdir"], "{EHA}_combined_contigs.fasta.gz.rev.2.bt2l"
            ),
        r1=os.path.join(
            config["workdir"], "reads/", "{PRB}/", "{EHI}_M_1.fq.gz"
            ),
        r2=os.path.join(
            config["workdir"], "reads/", "{PRB}/", "{EHI}_M_2.fq.gz"
            ),
        contigs=os.path.join(
            config["workdir"], "{EHA}_combined_contigs.fasta.gz"
            ),
    output:
        os.path.join(
            config["workdir"], "bams/", "{PRB}_{EHI}_{EHA}.bam"
            )
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 8
    resources:
        mem_gb=48,
        time=estimate_time_mapping,
    benchmark:
        os.path.join(config["logdir"] + "/multisplit_mapping_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/multisplit_mapping_log_{PRB}_{EHI}_{EHA}.log")
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