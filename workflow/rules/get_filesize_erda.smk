################################################################################
### Get file sizes from ERDA
rule filesize_from_ERDA:
    input:
        os.path.join(
            config["workdir"], 
            "ERDA_folder_created"
        )
    output:
        temp(
            os.path.join(
                config["workdir"],
                "{sample}_filesize.txt"
            )
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time='00:00:30'
    message:
        "Fetching filesize for {wildcards.sample} from ERDA"
    shell:
        """
        echo 'ls -l /EarthHologenomeInitiative/Data/RAW/*/{wildcards.sample}*_1.fq.gz' | sftp erda | sed '1d;' | awk '{{print $5}}' > {output}

        """