################################################################################
### Submit MAG information to ENA
rule register_mags_ena:
    input:
        sample_receipt=os.path.join(
            config["workdir"],
            "{EHI}_sample_receipt.tsv"
        )
    output:
        mag_receipt=os.path.join(
            config["workdir"],
            "{EHI}_mag_receipt.tsv"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time=00:05:00
    message:
        "Registering MAGs at the ENA for {wildcards.EHI}"
    shell:
        """

        ena-upload
        
        """