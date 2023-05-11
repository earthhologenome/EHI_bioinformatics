###############################################################################
## Upload tables to AirTable MAG database
rule upload_tables:
    input:
        count_table=os.path.join(
            config["workdir"], 
            "coverm/", 
            config["dmb"] + "_count_table.tsv"
        ),
        mapping_rates=os.path.join(
            config["workdir"], 
            "coverm/", 
            config["dmb"] + "_mapping_rate.tsv"
        ),
        combined=os.path.join(
            config["workdir"], 
            config["dmb"] + "_gtdbtk_combined_summary.tsv"
        )        
    output:
        os.path.join(
            config["workdir"],
            "tables_uploaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='00:10:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_tables_benchmark.tsv")    
    shell:
        """
        ## Upload MAGs and annotations to ERDA
        lftp sftp://erda -e "mirror -R {config[magdir]} /EarthHologenomeInitiative/Data/MAG/; bye"

        ## Add MAG mapping rates to airtable
        python {config[codedir]}/airtable/add_mag_mapping_rates_airtable.py --table={input.mapping_rates}

        ## Upload count table to airtable
        python {config[codedir]}/airtable/upload_count_table_airtable.py --table={input.count_table}

        ## Clean up

        ## Create output to end pipeline
        touch {output}
        """