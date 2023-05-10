################################################################################
### Estimate diversity and required sequencing effort using nonpareil
rule nonpareil:
    input:
        non_host_r1=os.path.join(
            config["workdir"],
            "{sample}_M_1.fq"
        ),
        non_host_r2=os.path.join(
            config["workdir"],
            "{sample}_M_2.fq"
        )
    output:
        npo=os.path.join(
            config["workdir"],
            "misc/{sample}.npo"
        ),
        npstats=os.path.join(
            config["workdir"],
            "misc/{sample}_np.tsv"
        ),
        non_host_r1c=temp(
            os.path.join(
                config["workdir"],
                "{sample}_M_1.fq.gz"
            )
        ),
        non_host_r2c=temp(
            os.path.join(
                config["workdir"],
                "{sample}_M_2.fq.gz"
            )
        )
    conda:
        f"{config['codedir']}/conda_envs/nonpareil.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=8,
        time=estimate_time_nonpareil
    benchmark:
        os.path.join(config["logdir"] + "/{sample}_nonpareil.benchmark.tsv")
    message:
        "Estimating microbial diversity using nonpareil"
    shell:
        """
        #IF statement to account for situations where there are not enough
        #microbial reads in a sample (e.g. high host% or non-metagenomic sample)
        #In this case, if R1 has > 150 Mbytes, run, else, skip:
        if [ $(( $(stat -c '%s' {input.non_host_r1}) / 1024 / 1024 )) -gt 150 ]
        then
        #Run nonpareil
        nonpareil \
            -s {input.non_host_r1} \
            -f fastq \
            -T kmer \
            -t {threads} \
            -b {wildcards.sample}

        mv {wildcards.sample}.* {config[workdir]}/misc/

        #Script to extract nonpareil values of interest
        Rscript {config[codedir]}/scripts/nonpareil_table.R {output.npo} {output.npstats}

        else
        #Create dummy file for snakemake to proceed
        touch {output.npo}
        echo -e "sample\tkappa\tC\tLR\tmodelRt\tLRstar\tdiversity\n{wildcards.sample}\t0\t0\t0\t0\t0\t0" > {output.npstats}

        fi


        #Compress reads
        pigz -p {threads} {input.non_host_r1}
        pigz -p {threads} {input.non_host_r2}

        """