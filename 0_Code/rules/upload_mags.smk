###############################################################################
## Upload MAGs to ERDA and update AirTable MAG database
rule upload_mags:
    input:
        expand(
            os.path.join(
                config["magdir"],
                "{MAG}_anno.tsv.gz"
            ),
            MAG=MAG
        )
    output:
            os.path.join(
                config["magdir"],
                "MAGs_uploaded"
                )
    shell:
        """
        ## Upload MAGs and annotations to ERDA
        lftp sftp://erda -e "mirror -R {config[magdir]} -o /EarthHologenomeInitiative/Data/MAG/; bye"

        ## Update AirTable MAG database
#        python {config[codedir]}/airtable/add_mags_airtable.py --report={output} --asb={config[abb]}

        ## Create output to end pipeline
        touch {output}
        """