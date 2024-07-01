################################################################################
### Log AirTable that the pipeline has started.
rule log_start:
    output:
        os.path.join(config["workdir"], "pipeline_started")
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time="00:05:00",
    message:
        "Logging AirTable that the run has started."
    shell:
        """
        python {config[codedir]}/airtable/log_ehs_start_airtable.py  --ehs={config[ehs]}

        touch {output}
        """