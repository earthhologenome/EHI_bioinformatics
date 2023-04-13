################################################################################
### Upload assembly BAMs to ERDA
rule upload_bam_erda:
    input:
        os.path.join(config["workdir"], "{PRB}/", "{EHI}/", "{EHI}_{EHA}.bam")
    output:
        os.path.join(config["workdir"], "{PRB}/", "{EHI}/", "{EHI}_{EHA}_uploaded")
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='00:30:00'
    benchmark:
        os.path.join(config["logdir"] + "upload_bam_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    message:
        "Uploading {wildcards.EHA} BAM to ERDA"
    shell:
        """
        #Upload BAM to ERDA for storage
        lftp sftp://erda -e "put {input} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
        touch {output}
        """