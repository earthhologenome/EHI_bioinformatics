################################################################################
### Log AirTable that the pipeline has finished.
### Also moves MAGs for the 3_mag_annotation.snakefile
rule log_finish:
    input:
        expand(
            os.path.join(
                "/projects/ehi/data/REP/",
                "{combo[0]}_{combo[1]}_{combo[2]}_final_stats.tsv",
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
        python {config[codedir]}/airtable/log_asb_done_airtable.py  --code={config[abb]}

        # Move MAGs for the next pipeline
        mkdir -p {config[magdir]}

        for mag in {config[workdir]}/*_refinement/metawrap_50_10_bins/*.fa.gz;
            do mv $mag {config[magdir]}/;
        done

        # Upload MAGs and annotations to ERDA
        lftp sftp://erda -e "mirror -R {config[magdir]} /EarthHologenomeInitiative/Data/MAG/; bye"

        # Clean up
        rm -r {config[magdir]}/*

        # Clean up
        rm -r {config[workdir]}/*

        touch {output}
        """