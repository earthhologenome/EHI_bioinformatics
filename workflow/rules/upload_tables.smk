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
            config["dmb"] + ".tree"
        ),
        counts=os.path.join(
            config["workdir"],
            config["dmb"] + "_counts.tsv"
        ),
        coverage=os.path.join(
            config["workdir"],
            config["dmb"] + "_coverage.tsv"
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

        python {config[codedir]}/airtable/get_sample_metadata_airtable.py --samples=read_input.tsv --dmb={config[dmb]}

        ## Remove trailing '.fa'
        sed -i'' 's/\.fa//g' {config[dmb]}_mag_info.tsv

        ## Upload other files to AirTable (count table, tree)
        gzip -k {input.count_table}
        gzip -k {config[dmb]}_metadata.tsv
        gzip -k {input.tree}
        gzip -k {input.pruned_tree}
        gzip -k {config[dmb]}_mag_info.tsv
        gzip -k {input.combined}
        gzip -k {input.counts}
        gzip -k {input.coverage}

        lftp sftp://erda -e "put {input.count_table}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.tree}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.pruned_tree}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.combined}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {config[dmb]}_mag_info.tsv.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {config[dmb]}_metadata.tsv.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.counts}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.coverage}.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {config[workdir]}/mag_catalogue/{config[dmb]}_mags.fasta.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 5
        lftp sftp://erda -e "mirror -R {config[workdir]}/bams/ /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"
        sleep 60

        ## Log # of dereplicated MAGs to airtable
        ls -l {config[workdir]}/drep/dereplicated_genomes/*.fa.gz | wc -l > {config[dmb]}_dereplicated_mags.tsv

        ## Log dereplicated MAG to airtable
        python {config[codedir]}/airtable/log_dmb_derep_names_airtable.py --mags={config[workdir]}/dereplicated_mags.tsv

        ## Log AirTable that the run is finished
        python {config[codedir]}/airtable/log_dmb_done_airtable.py --code={config[dmb]}

        ## Clean up
#        rm -r {config[workdir]}/*

        ## Create dummy files to end pipeline
        touch {output}
        """