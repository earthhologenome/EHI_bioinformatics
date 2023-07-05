###############################################################################
## Upload DRAM annotations to ERDA and update AirTable MAG database
rule upload_mags:
    input:
        mags=expand(
            os.path.join(config["magdir"], "{ehm}_anno.tsv.gz"),
            ehm=[combo[0] for combo in valid_combinations]
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
        time='04:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_mag_benchmark.tsv")    
    shell:
        """
        #Setup batch file for uploading MAGs from erda:
        for mag in {input.mags};
            do echo "put" >> {config[workdir]}/put.tsv && echo $(basename $mag) >> {config[workdir]}/up_mag.tsv && echo "erda:EarthHologenomeInitiative/Data/ANN/" > ann.tsv;
        done

        paste {config[workdir]}/get.tsv {config[workdir]}/up_mag.tsv {config[workdir]}/ann.tsv -d '' > {config[workdir]}/upload_batchfile.txt

        #Execute batch file to upload the suckers
        cd {config[magdir]}
        sftp -b {config[workdir]}/upload_batchfile.txt erda

        ## Clean up
#        rm -r {config[magdir]}/*

        ## Log job is done on AirTable
        python {config[codedir]}/airtable/log_ann_done_airtable.py --dmb={config[dmb]}

        ## Create output to end pipeline
        touch {output}
        """