################################################################################
### Calculate % of each sample's reads mapping to host genome/s (also upload PPR reads to ERDA)
rule coverM:
    input:  
        npo=os.path.join(
            config["workdir"],
            "misc/{sample}.npo"
        ),
        bam=os.path.join(
            config["workdir"],
            "tmp/{sample}.bam"
        )
    output:
        os.path.join(
            config["workdir"],
            "misc/{sample}_coverM_mapped_host.tsv"
        )
    conda:
        f"{config['codedir']}/conda_envs/coverm.yaml"
    threads:
        2
    resources:
        load=1,
        mem_gb=8,
        time='00:10:00'
    benchmark:
        os.path.join(config["logdir"] + "/{sample}_coverM.benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/{sample}_coverM.log")
    message:
        "Calculating percentage of reads mapped to host genome/s using coverM"
    shell:
        """
        #Calculate % mapping to host using coverM
        coverm genome \
            -b {input.bam} \
            -s _ \
            -m relative_abundance count \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}

        #Remove empty nonpareil (.npo) files to streamline future plotting
        if [ $(stat -c '%s' {input.npo}) -lt 1 ]
        then
        rm {input.npo}
        fi

        """