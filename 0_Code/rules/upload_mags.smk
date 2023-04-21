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
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='01:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_mag_benchmark.tsv")    
    shell:
        """
        ## Upload MAGs and annotations to ERDA
        lftp sftp://erda -e "mirror -R {config[magdir]} /EarthHologenomeInitiative/Data/MAG/; bye"

        ## Clean up
        rm -r {config[magdir]}/*

        ## Create output to end pipeline
        touch {output}
        """