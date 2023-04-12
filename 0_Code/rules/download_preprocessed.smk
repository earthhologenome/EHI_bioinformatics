################################################################################
### Fetch preprocessed reads from ERDA
rule download_preprocessed_from_ERDA:
    output:
        r1=os.path.join(config["workdir"], "{PRB}", "{EHI}_M_1.fq.gz"),
        r2=os.path.join(config["workdir"], "{PRB}", "{EHI}_M_2.fq.gz"),
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time="00:15:00",
    message:
        "Fetching metagenomics reads for {wildcards.EHI} from ERDA"
    shell:
        """
        lftp sftp://erda -e "mirror --include-glob='{wildcards.PRB}/{wildcards.EHI}*.fq.gz' /EarthHologenomeInitiative/Data/PPR/ {{config['workdir']}}/; bye"
        """