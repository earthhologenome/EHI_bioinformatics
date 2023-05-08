################################################################################
### Calculate % of each sample's reads mapping to host genome/s (also upload PPR reads to ERDA)
rule upload_to_ERDA:
    input:
        pipe=os.path.join(
            config["workdir"],
            "misc/{sample}_pipe.tsv.gz"
        ),      
        non_host_r1=os.path.join(
            config["workdir"],
            "{sample}_M_1.fq.gz"
        ),
        non_host_r2=os.path.join(
            config["workdir"],
            "{sample}_M_2.fq.gz"
        ),
        host_bam=os.path.join(
            config["workdir"],
            "{sample}_G.bam"
        ),
        coverm=os.path.join(
            config["workdir"],
            "misc/{sample}_coverM_mapped_host.tsv"
        ),      
        file_size=os.path.join(
            config["workdir"],
            "{sample}_filesize.txt"
        )
    output:
        os.path.join(
            config["workdir"],
            "misc/{sample}_uploaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=16,
        time=estimate_time_download
    log:
        os.path.join(config["logdir"] + "/{sample}_upload.log")
    message:
        "Uploading reads and BAM to ERDA"
    shell:
        """
        #Upload preprocessed reads to ERDA for storage
        lftp sftp://erda -e "put {input.non_host_r1} -o /EarthHologenomeInitiative/Data/PPR/{config[prb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.non_host_r2} -o /EarthHologenomeInitiative/Data/PPR/{config[prb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.host_bam} -o /EarthHologenomeInitiative/Data/PPR/{config[prb]}/; bye"
        
        touch {output}
        """