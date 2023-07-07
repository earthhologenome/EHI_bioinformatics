###############################################################################
## Upload DRAM annotations to ERDA and update AirTable MAG database
rule upload_mags:
    input:
        mags=expand(
            os.path.join(
                config["magdir"], "{MAG}_anno.tsv.gz"
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
        time='04:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_mag_benchmark.tsv")    
    shell:
        """
        ##Rename files from EHA -> EHM
        sed -s '1d;' dereped_mags.csv | tr ',' '\t' > ehm_eha_mapping.tsv
        #fix issue with separators
        dos2unix ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha"_anno.tsv.gz {config[magdir]}/"$ehm"_anno.tsv.gz && echo {config[magdir]}/"$ehm"_anno.tsv.gz >> anno_mag.tsv; 
        done < ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha"_kegg.tsv.gz {config[magdir]}/"$ehm"_kegg.tsv.gz && echo {config[magdir]}/"$ehm"_kegg.tsv.gz >> kegg_mag.tsv; 
        done < ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha".gbk.gz {config[magdir]}/"$ehm".gbk.gz && echo {config[magdir]}/"$ehm".gbk.gz >> gbk_mag.tsv; 
        done < ehm_eha_mapping.tsv


        #Setup batch file for uploading MAGs from erda:
        for mag in {config[magdir]}/EHM*_anno.tsv.gz;
            do echo "put" >> put.tsv && echo "EarthHologenomeInitiative/Data/ANN/" >> ann.tsv;
        done

        cat anno_mag.tsv kegg_mag.tsv gbk_mag.tsv > upload_filenames.tsv
        cat put.tsv put.tsv put.tsv > upload_put.tsv
        cat ann.tsv ann.tsv ann.tsv > upload_ann.tsv

        paste upload_put.tsv upload_filenames.tsv upload_ann.tsv > batchfile.txt

        #Execute batch file to upload the suckers
        sftp -b batchfile.txt erda

        ## Clean up
        rm -r {config[magdir]}/*

        ## Log job is done on AirTable
        python {config[codedir]}/airtable/log_ann_done_airtable.py --code={config[dmb]}

        ## Create output to end pipeline
        touch {output}
        """