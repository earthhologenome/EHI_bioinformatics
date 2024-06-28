################################################################################
### Upload MAGs to ENA using Aspera
rule upload_mags_aspera:
    input:
        mag_receipt=os.path.join(
            config["workdir"],
            "{EHI}_sample_receipt.tsv"
        )
    output:
        mags_uploaded=os.path.join(
            config["workdir"],
            "{EHI}_mags_uploaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time=00:05:00
    message:
        "Uploading MAGs to ENA for {wildcards.EHI}"
    shell:
        """
        aspera

        touch {output.mags_uploaded}
        
        """