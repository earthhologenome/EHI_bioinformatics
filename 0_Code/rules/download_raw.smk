################################################################################
### Fetch raw data from ERDA
rule download_from_ERDA:
    input:
        ready=os.path.join(
            config["workdir"], 
            "ERDA_folder_created"
        ),
        filesize=os.path.join(
            config["workdir"],
            "{sample}_filesize.txt"
        )
    output:
        r1o=temp(
            os.path.join(
                config["workdir"],
                "{sample}_raw_1.fq.gz"
            )
        ),
        r2o=temp(
            os.path.join(
                config["workdir"],
                "{sample}_raw_2.fq.gz"
            )
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time=estimate_time_download
    message:
        "Fetching {wildcards.sample} from ERDA"
    shell:
        """
        lftp sftp://erda -e "mirror --include-glob='{wildcards.sample}*.fq.gz' /EarthHologenomeInitiative/Data/RAW/ {config[workdir]}; bye"
        mv {config[workdir]}/SEB*/{wildcards.sample}*_1.fq.gz {output.r1o}
        mv {config[workdir]}/SEB*/{wildcards.sample}*_2.fq.gz {output.r2o}
        """