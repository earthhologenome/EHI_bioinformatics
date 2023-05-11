################################################################################
### Create EHA folder on ERDA
rule create_ASB_folder:
    output:
        os.path.join(
            config["workdir"], 
            "ERDA_folder_created"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time="00:05:00",
    message:
        "Creating assembly batch folder on ERDA"
    shell:
        """
        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/ASB/{config[abb]} ; bye"
        touch {output}

        #Also, log the AirTable that the ASB is running!
        python {config[codedir]}/airtable/log_asb_start_airtable.py --code={config[abb]}
        """