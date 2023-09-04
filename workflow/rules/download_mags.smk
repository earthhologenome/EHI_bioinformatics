################################################################################
### Fetch MAGs from ERDA
rule download_mags:
    input:
        os.path.join(
            config["workdir"], 
            "ERDA_folder_created"
        )
    output:
        mags=expand(
            os.path.join(
                config["magdir"], 
                "{mag}.gz"
            ), mag = MAG
        ),
        downloaded=os.path.join(
            config["magdir"],
            "mags_downloaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time=estimate_time_download_mags
    benchmark:
        os.path.join(config["logdir"] + "/download_mags_benchmark.tsv")
    message:
        "Fetching MAGs from ERDA"
    shell:
        """
        rm -f {config[workdir]}/batchfile.txt
        rm -f {config[workdir]}/mag.tsv
        rm -f {config[workdir]}/get.tsv

        #Setup batch file for downloading MAGs from erda:
        for mag in {output.mags};
            do echo "get EarthHologenomeInitiative/Data/MAG/*/" >> {config[workdir]}/get.tsv && echo $(basename $mag) >> {config[workdir]}/mag.tsv;
        done

        paste {config[workdir]}/get.tsv {config[workdir]}/mag.tsv -d '' > {config[workdir]}/batchfile.txt

        #Execute batch file to pull the suckers
        cd {config[magdir]}
        sftp -b {config[workdir]}/batchfile.txt erda

        #Indicate files are downloaded
        touch {output.downloaded}
        """