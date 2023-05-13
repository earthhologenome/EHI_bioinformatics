################################################################################
### Map reads to MAG catalogue
rule mag_mapping:
    input:
        contigs=os.path.join(
            config["workdir"],
            "mag_catalogue/",
            config["dmb"] + "_mags.fasta.gz"
            ),
        index=os.path.join(
            config["workdir"],
            "mag_catalogue/",
            config["dmb"] + "_mags.fasta.gz.rev.2.bt2l"
            ),
        r1=os.path.join(
            config["workdir"], 
            "reads/", 
            "{PRB}/", 
            "{EHI}_M_1.fq.gz"
            ),
        r2=os.path.join(
            config["workdir"], 
            "reads/", 
            "{PRB}/", 
            "{EHI}_M_2.fq.gz"
            )
    output:
        os.path.join(
            config["workdir"], 
            "bams/", 
            "{PRB}_{EHI}_" + config["dmb"] + ".bam"
            )
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 8
    resources:
        mem_gb=48,
        time="04:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/mag_mapping_benchmark_{PRB}_{EHI}.tsv")
    log:
        os.path.join(config["logdir"] + "/mag_mapping_log_{PRB}_{EHI}.log")
    message:
        "Mapping {wildcards.EHI} to MAG catalogue using Bowtie2"
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