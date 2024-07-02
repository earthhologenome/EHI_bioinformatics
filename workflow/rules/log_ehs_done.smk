################################################################################
### Log AirTable that the pipeline has finished.
rule log_finish:
    input:
        expand(
            os.path.join(
            config["workdir"],
            "{EHI}_accessions_uploaded"
            ),
            EHI=EHI
        )
    output:
        os.path.join(
            config["workdir"], 
            "pipeline_finished"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time="00:05:00",
    message:
        "Logging AirTable that the pipeline has been completed."
    shell:
        """
        # Log on the AirTable that the pipeline has finished:
        python {config[codedir]}/airtable/log_ehs_done_airtable.py  --ehs={config[ehs]}

        touch {output}
        """