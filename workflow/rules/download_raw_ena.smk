################################################################################
### Fetch raw data from ERDA
rule download_raw:
    output:
        r1=temp(
            os.path.join(
                config["workdir"],
                "{EHI}_raw_1.fq.gz"
            )
        ),
        r2=temp(
            os.path.join(
                config["workdir"],
                "{EHI}_raw_2.fq.gz"
            )
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time=estimate_time_download
    message:
        "Fetching {wildcards.EHI} from ERDA"
    shell:
        """
        wget --no-verbose `grep '{wildcards.EHI}' ehi_numbers.tsv | cut -f3 | sed 's/"//g'`
        mv {wildcards.EHI}*.fq.gz {output.r1}

        wget --no-verbose `grep '{wildcards.EHI}' ehi_numbers.tsv | cut -f4 | sed 's/"//g'`
        mv {wildcards.EHI}*.fq.gz {output.r2}
        """