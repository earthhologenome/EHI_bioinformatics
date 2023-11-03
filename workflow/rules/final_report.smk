###############################################################################
## Create and upload final report
rule final_report:
    input:
        os.path.join(
            config["magdir"],
            "MAGs_uploaded"
        )
    output:
        os.path.join(
            config["magdir"],
            "report_uploaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    params:
        stats_dir=os.path.join(config["workdir"] + "/ANN_stats_dir")
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='01:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/final_report_benchmark.tsv")    
    shell:
        """
        #Activate conda environment
        source activate /projects/ehi/data/0_Environments/conda/final_report

        #Pull required input files for R script:
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}_counts.tsv.gz
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}_coverage.tsv.gz
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}_metadata.tsv.gz
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}_mag_info.tsv.gz
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}.tree.gz
        wget https://sid.erda.dk/share_redirect/BaMZodj9sA/DMB/{config[dmb]}/{config[dmb]}_merged_kegg.tsv.gz


        #Execute R script
        cp {config[codedir]}/scripts/final_report_render.r . 
        cp {config[codedir]}/scripts/final_report.r . 
        Rscript final_report_render.r -s final_report.r -b {config[dmb]} -p {config[dmb]}.pdf

        #Unload conda environment
        conda deactivate /projects/ehi/data/0_Environments/conda/final_report

        #Upload pdf to ERDA
        lftp sftp://erda -e "put {config[dmb]}.pdf -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"

        ## Clean up
        rm -r {config[magdir]}/*

        ## Log job is done on AirTable
        python {config[codedir]}/airtable/log_ann_done_airtable.py --code={config[dmb]}

        ## Create output to end pipeline
        touch {output}

        """