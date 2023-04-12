################################################################################
### Log AirTable that the pipeline has finished.
rule log_finish:
    input:
        expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}/",
                "{combo[2]}_final_stats.tsv",
            ),
            combo=valid_combinations,
        ),
    output:
        os.path.join(config["workdir"], "{abb}_pipeline_finished")
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time="00:05:00",
    message:
        "Logging AirTable that the run has been completed."
    shell:
        """
        # Log on the AirTable that the pipeline has finished:
        python {{config['codedir']}}/airtable/log_asb_done_airtable.py  --code={{config['abb']}}
        """