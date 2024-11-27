################################################################################
### Fetch dereplicated MAGs from ERDA
rule download_mags:
    output:
        mags=expand(
            os.path.join(
                config["magdir"], "{MAG}.gz"
                ),
                MAG=MAG
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
        time="08:00:00"
    benchmark:
        os.path.join(config["logdir"] + "/download_mags_benchmark.tsv")
    message:
        "Fetching MAGs from ERDA"
    shell:
        """
        #Log airtable that pipeline is running
        python {config[codedir]}/airtable/log_ann_start_airtable.py --code={config[dmb]}
        
        rm -rf {config[workdir]}/*

        #Strip characters from derep mag list
        sed '1d;' /projects/ehi/data/RUN/{config[dmb]}/dereped_mags.csv > {config[workdir]}/dereped_mags.csv

        sed 's/\[//g' dereped_mags.csv | sed 's/\]//g' | sed "s/'//g" | tr ',' '\t' > dereped_mags_clean.csv

        while read ehm eha abb; 
            do echo -e "get EarthHologenomeInitiative/Data/MAG/""$abb" >> {config[workdir]}/get.tsv
        done < dereped_mags_clean.csv

        while read ehm eha abb; 
            do echo -e "$eha"".gz" >> {config[workdir]}/eha.tsv
        done < dereped_mags_clean.csv

        dos2unix {config[workdir]}/eha.tsv && dos2unix {config[workdir]}/get.tsv 

        paste {config[workdir]}/get.tsv {config[workdir]}/eha.tsv -d '/' > {config[workdir]}/batchfile.txt

        #Execute batch file to pull the suckers
        mkdir -p {config[magdir]}
        cd {config[magdir]}
        sftp -b {config[workdir]}/batchfile.txt erda

        #Indicate files are downloaded
        touch {output.downloaded}
        """