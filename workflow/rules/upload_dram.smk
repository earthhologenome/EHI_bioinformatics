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
    params:
        stats_dir=os.path.join(config["workdir"] + "/ANN_stats_dir")
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='04:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/upload_mag_benchmark.tsv")    
    shell:
        """
        rm -rf {params.stats_dir}
        mkdir -p {params.stats_dir}

        ##Combine kegg outputs to a single table per DMB
        for i in {config[magdir]}/*_kegg.tsv.gz;
            do cat $i >> {params.stats_dir}/merged_kegg.tsv.gz;
        done

        gunzip {params.stats_dir}/merged_kegg.tsv.gz

        grep -v 'genome' {params.stats_dir}/merged_kegg.tsv > {params.stats_dir}/merged_kegg_body.tsv
        head -1 {params.stats_dir}/merged_kegg.tsv > {params.stats_dir}/merged_kegg_head.tsv

        cat {params.stats_dir}/merged_kegg_head.tsv {params.stats_dir}/merged_kegg_body.tsv > {params.stats_dir}/{config[dmb]}_merged_kegg.tsv
        gzip {params.stats_dir}/{config[dmb]}_merged_kegg.tsv

        lftp sftp://erda -e "put {params.stats_dir}/{config[dmb]}_merged_kegg.tsv.gz -o /EarthHologenomeInitiative/Data/DMB/{config[dmb]}/; bye"

        ##Rename files from EHA -> EHM
        sed -s '1d;' dereped_mags.csv | tr ',' '\t' > {params.stats_dir}/ehm_eha_mapping.tsv
        #fix issue with separators
        module load dos2unix/7.4.2
        dos2unix {params.stats_dir}/ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha"_anno.tsv.gz {config[magdir]}/"$ehm"_anno.tsv.gz && echo {config[magdir]}/"$ehm"_anno.tsv.gz >> {params.stats_dir}/anno_mag.tsv; 
        done < {params.stats_dir}/ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha"_kegg.tsv.gz {config[magdir]}/"$ehm"_kegg.tsv.gz && echo {config[magdir]}/"$ehm"_kegg.tsv.gz >> {params.stats_dir}/kegg_mag.tsv; 
        done < {params.stats_dir}/ehm_eha_mapping.tsv

        while read ehm eha; 
            do cp {config[magdir]}/"$eha".gbk.gz {config[magdir]}/"$ehm".gbk.gz && echo {config[magdir]}/"$ehm".gbk.gz >> {params.stats_dir}/gbk_mag.tsv; 
        done < {params.stats_dir}/ehm_eha_mapping.tsv   

        #Setup batch file for uploading MAGs from erda:
        for mag in {config[magdir]}/EHM*_anno.tsv.gz;
            do echo "put" >> {params.stats_dir}/put.tsv && echo "EarthHologenomeInitiative/Data/ANN/" >> {params.stats_dir}/ann.tsv;
        done

        cat {params.stats_dir}/anno_mag.tsv {params.stats_dir}/kegg_mag.tsv {params.stats_dir}/gbk_mag.tsv > {params.stats_dir}/upload_filenames.tsv
        cat {params.stats_dir}/put.tsv {params.stats_dir}/put.tsv {params.stats_dir}/put.tsv > {params.stats_dir}/upload_put.tsv
        cat {params.stats_dir}/ann.tsv {params.stats_dir}/ann.tsv {params.stats_dir}/ann.tsv > {params.stats_dir}/upload_ann.tsv

        paste {params.stats_dir}/upload_put.tsv {params.stats_dir}/upload_filenames.tsv {params.stats_dir}/upload_ann.tsv > {params.stats_dir}/batchfile.txt

        #Execute batch file to upload the suckers
        sftp -b {params.stats_dir}/batchfile.txt erda

        ## Log DRAM results in AirTable
        for i in {config[magdir]}/*annotate;
            do echo $(basename ${{i/_annotate/}}) >> {params.stats_dir}/mag_names_at.tsv && cat $i/cazy_hits.tsv >> {params.stats_dir}/cazy_hits_at.tsv && cat $i/pfam_hits.tsv >> {params.stats_dir}/pfam_hits_at.tsv && cat $i/kegg_hits.tsv >> {params.stats_dir}/kegg_hits_at.tsv && cat $i/unannotated.tsv >> {params.stats_dir}/unannotated_at.tsv && cat $i/num_genes.tsv >> {params.stats_dir}/num_genes_at.tsv;
        done

        echo -e "mag_name\tnumber_genes\tcazy_hits\tpfam_hits\tkegg_hits\tnumber_unannotated_genes" > {params.stats_dir}/at_headers.tsv
        paste {params.stats_dir}/mag_names_at.tsv {params.stats_dir}/num_genes_at.tsv {params.stats_dir}/cazy_hits_at.tsv {params.stats_dir}/pfam_hits_at.tsv {params.stats_dir}/kegg_hits_at.tsv {params.stats_dir}/unannotated_at.tsv > {params.stats_dir}/dram_at.tsv
        cat {params.stats_dir}/at_headers.tsv {params.stats_dir}/dram_at.tsv > {params.stats_dir}/at_dram.tsv

        python {config[codedir]}/airtable/add_mag_dram_results_airtable.py --table={params.stats_dir}/at_dram.tsv

        ## Create output to trigger rule end
        touch {output}

        """