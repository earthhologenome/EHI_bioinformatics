################################################################################
### Calculate the number of reads that mapped to the mag catalogue
rule coverM_mag:
    input:
        bam=expand(
            os.path.join(
                config["workdir"], 
                "bams/", 
                "{combo[0]}_{combo[1]}_" + config["dmb"] + ".bam"
            ), combo=valid_combinations
        )
    output:
        count_table=os.path.join(
            config["workdir"], 
            "coverm/", 
            config["dmb"] + "_count_table.tsv"
        ),
        mapping_rate=os.path.join(
            config["workdir"], 
            "coverm/", 
            config["dmb"] + "_mapping_rate.tsv"
        )
    conda:
        f"{config['codedir']}/conda_envs/coverm.yaml"
    threads: 4
    resources:
        load=1,
        mem_gb=64,
        time="00:15:00",
    benchmark:
        os.path.join(config["logdir"] + "/coverm_benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/coverm_log.log")
    message:
        "Calculating mag mapping rate with CoverM"
    shell:
        """
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_fraction length \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.count_table}

        #relative abundance for report
        coverm genome \
            -b {input} \
            -s ^ \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.mapping_rate}

        """