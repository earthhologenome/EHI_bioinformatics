###############################################################################
## Upload MAGs to ERDA and update AirTable MAG database
rule upload_mags:
    input:
        expand(
            os.path.join(
                config["workdir"],
                "{PRB}/",
                "{EHI}/",
                "{EHA}/",
                "DRAM/",
                "{MAG}_anno.tsv.gz"
                ),
        )
    output:
            os.path.join(
                config["workdir"],
                "{PRB}/",
                "{EHI}/",
                "{EHA}/",
                "DRAM/",
                "MAGs_uploaded"
                )
    shell:
        """
        ## Upload MAGs and annotations to ERDA
        lftp sftp://erda -e "mirror -R {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}/DRAM/ -o /EarthHologenomeInitiative/Data/MAG/{wildcards.EHA}/; bye"

        ## Update AirTable MAG database
        python {config[codedir]}/airtable/add_mags_airtable.py --report={output} --asb={config[abb]}

        touch {output}
        """