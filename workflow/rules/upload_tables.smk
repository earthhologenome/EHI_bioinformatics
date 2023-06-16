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
        ),
        tree=os.path.join(
            config["workdir"],
            config["dmb"] + "_gtdbtk.bac120.classify.tree"
        ),
        pruned_tree=os.path.join(
            config["workdir"],
            config["dmb"] + "_pruned.tree"
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
        time='04:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_tables_benchmark.tsv")    
    shell:
        """
        ## Add MAG mapping rates to airtable
        head -2 {input.mapping_rates} | cut -f2- | sed 's/ Relative Abundance (%)//g' > unmapped.tsv

        # Print the header
        echo -e "PR_batch_static\tEHI_sample_static\tDM_batch_static\tMAG_mapping_percentage" > mapping_header.tsv

        # Run script to transpose table
        bash {config[codedir]}/scripts/transpose_table.sh

        cat mapping_header.tsv longer.tsv > mapping_rates.tsv
        
        python {config[codedir]}/airtable/add_mag_mapping_rates_airtable.py --report=mapping_rates.tsv

        python {config[codedir]}/airtable/get_mag_info_airtable.py --dmb={config[dmb]}

        ## Upload other files to AirTable (count table, tree)
        gzip {input.count_table}
        gzip {input.tree}
        gzip {input.pruned_tree}
        gzip {config[dmb]}_mag_info.tsv
        gzip {input.combined}


        lftp sftp://erda -e "put {input.count_table}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.tree}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.pruned_tree}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.combined}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {config[dmb]}__mag_info.tsv -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {config[workdir]}/bams/*.bam -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 60

        ## Log # of dereplicated MAGs to airtable
        ls -l {config[workdir]}/drep/dereplicated_genomes/*.fa.gz | wc -l > dereplicated_mags.tsv


        ## Log AirTable that the run is finished
        python {config[codedir]}/airtable/log_dmb_done_airtable.py --code={config[dmb]}

        ## Clean up
        rm -r {config[workdir]}/*

        ## Create output to end pipeline
        touch {output}
        """