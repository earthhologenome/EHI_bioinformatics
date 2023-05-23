################################################################################
### Upload assembly BAMs to ERDA
rule upload_bam_erda:
    input:
        bam=os.path.join(
            config["workdir"], 
            "bams/" 
            "{PRB}_{EHI}_{EHA}.bam"
        ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
        )
    output:
        os.path.join(config["workdir"], "{PRB}_{EHI}_{EHA}_uploaded")
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time=estimate_time_upload_bam
    benchmark:
        os.path.join(config["logdir"] + "/upload_bam_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    message:
        "Uploading {wildcards.EHA} BAM to ERDA"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt 40 ]
        then
            touch {output}

        else

            #Upload BAM to ERDA for storage
            lftp sftp://erda -e "put {input.bam} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
            touch {output}
        fi
        """