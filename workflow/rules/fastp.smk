################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i=os.path.join(
            config["workdir"],
            "{sample}_raw_1.fq.gz"
        ),
        r2i=os.path.join(
            config["workdir"],
            "{sample}_raw_2.fq.gz"
        )
    output:
        r1o=temp(
            os.path.join(
                config["workdir"],
                "{sample}_trimmed_1.fq.gz"
            )
        ),
        r2o=temp(
            os.path.join(
                config["workdir"],
                "{sample}_trimmed_2.fq.gz"
            )
        ),
        fastp_html=os.path.join(
            config["workdir"],
            "misc/{sample}.html"
        ),
        fastp_json=os.path.join(
            config["workdir"],
            "misc/{sample}.json"
        )
    params:
        adapter1=expand("{adapter1}", adapter1=config['adapter1']),
        adapter2=expand("{adapter2}", adapter2=config['adapter2'])
    conda:
        f"{config['codedir']}/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=10,
        time=estimate_time_fastp
    benchmark:
        os.path.join(config["logdir"] + "/{sample}_fastp.benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/{sample}_fastp.log")
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
            --low_complexity_filter \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence {params.adapter1} \
            --adapter_sequence_r2 {params.adapter2} \
        &> {log}
        """