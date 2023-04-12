################################################################################
### Map reads to assemblies
rule assembly_mapping:
    input:
        index=os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_contigs.fasta.rev.2.bt2l"
        ),
        r1=os.path.join(config["workdir"], "{PRB}", "{EHI}_M_1.fq.gz"),
        r2=os.path.join(config["workdir"], "{PRB}", "{EHI}_M_2.fq.gz"),
        contigs=os.path.join(config["workdir"], "{PRB}" "{EHI}" "{EHA}_contigs.fasta"),
    output:
        os.path.join(config["workdir"], "{PRB}", "{EHI}", "{EHI}_{EHA}.bam"),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=48,
        time="05:00:00",
    benchmark:
        "{{config['logdir']}}/assembly_mapping_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/assembly_mapping_log_{PRB}_{EHI}_{EHA}.log",
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