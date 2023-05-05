################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i=os.path.join(
            config["workdir"],
            "{sample}_trimmed_1.fq.gz"
        ),
        r2i=os.path.join(
            config["workdir"],
            "{sample}_trimmed_2.fq.gz"
        ),
        bt2_index=os.path.join(
            config["workdir"],
            config["hostgenome"],
            config["hostgenome"] + "_RN.fna.gz.rev.2.bt2l",
        ),
        catted_ref=os.path.join(
            config["workdir"],
            config["hostgenome"],
            config["hostgenome"] + "_RN.fna.gz"
        )
    output:
        all_bam=temp(os.path.join(
            config["workdir"],
            "tmp/{sample}.bam"
            )
        ),
        host_bam=temp(os.path.join(
            config["workdir"],
            "{sample}_G.bam"
            )
        ),
        non_host_r1=os.path.join(
            config["workdir"],
            "{sample}_M_1.fq"
        ),
        non_host_r2=os.path.join(
            config["workdir"],
            "{sample}_M_2.fq"
        )
    conda:
        f"{config['codedir']}/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=24,
        time=estimate_time_mapping
    benchmark:
        os.path.join(config["logdir"] + "/{sample}_mapping.benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/{sample}_mapping.log")
    message:
        "Mapping {wildcards.sample} reads to host genomes"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.catted_ref} \
            -1 {input.r1i} \
            -2 {input.r2i} \
        | samtools view -b -@ {threads} - | samtools sort -@ {threads} -o {output.all_bam} - &&

        # Extract non-host reads (note we're not compressing for nonpareil)
        samtools view -b -f12 -@ {threads} {output.all_bam} \
        | samtools fastq -@ {threads} -1 {output.non_host_r1} -2 {output.non_host_r2} - &&

        # Send host reads to BAM
        samtools view -b -F12 -@ {threads} {output.all_bam} \
        | samtools sort -@ {threads} -o {output.host_bam} -
        """