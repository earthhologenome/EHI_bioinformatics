###############################################################################
## Upload DRAM annotations to ERDA and update AirTable MAG database
rule upload_mags:
    input:
        mags=expand(
            os.path.join(
                config["magdir"], "{MAG}_anno.tsv.gz"
                ),
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

        while read ehm eha; 
            do mv ${{eha/.fa/_anno.tsv.gz}} "$ehm"_anno.tsv.gz && echo "$ehm"_anno.tsv.gz >> anno_mag.tsv; 
        done < ehm_eha_mapping.tsv

        while read ehm eha; 
            do mv ${{eha/.fa/_kegg.tsv.gz}} "$ehm"_kegg.tsv.gz && echo "$ehm"_kegg.tsv.gz >> kegg_mag.tsv; 
        done < ehm_eha_mapping.tsv

        while read ehm eha; 
            do mv ${{eha/.fa/.gbk.gz}} "$ehm".gbk.gz && echo "$ehm"_gbk.gz >> gbk_mag.tsv; 
        done < ehm_eha_mapping.tsv


        #Setup batch file for uploading MAGs from erda:
        for mag in *_anno.tsv.gz;
            do echo "put" >> put.tsv && echo "erda:EarthHologenomeInitiative/Data/ANN/" > ann.tsv;
        done

        cat anno_mag.tsv kegg_mag.tsv gbk_mag.tsv > upload_filenames.tsv
        cat put.tsv put.tsv put.tsv > upload_put.tsv
        cat ann.tsv ann.tsv ann.tsv > upload_ann.tsv

        paste upload_filenames.tsv upload_put.tsv upload_ann.tsv > batchfile.txt

        #Execute batch file to upload the suckers
        sftp -b batchfile.txt erda

        ## Clean up
#        rm -r {config[magdir]}/*

        ## Log job is done on AirTable
        python {config[codedir]}/airtable/log_ann_done_airtable.py --dmb={config[dmb]}

        ## Create output to end pipeline
        touch {output}
        """