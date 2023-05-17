###############################################################################
## Generate output summary for coassemblies, update airtable
rule coassembly_summary:
    input:
        mw_stats=os.path.join(
            config["workdir"],
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        coverm=expand(os.path.join(
            config["workdir"], 
            "coverm/", 
            "{combo[0]}_{combo[1]}_{combo[2]}_assembly_coverM.txt"
                ), combo=valid_combinations
            ),
        tarball=expand(os.path.join(
            config["workdir"], 
            "coverm/", 
            "{combo[0]}_{combo[1]}_{combo[2]}_coverm.tar.gz"
                ), combo=valid_combinations
            ),
        contigs=os.path.join(
            config["workdir"], 
            "{EHA}_assembly/", 
            "{EHA}_contigs.fasta"
            ),
        gtdb=os.path.join(
            config["workdir"], 
            "{EHA}/",
            "{EHA}_gtdbtk_combined_summary.tsv"
            )
    output:
        stats=os.path.join(
            "/projects/ehi/data/REP/",
            "{EHA}_final_stats.tsv",
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    params:
        quast=directory(os.path.join(
            config["workdir"], 
            "{EHA}_QUAST")
            ),
        stats_dir=directory(os.path.join(
            config["workdir"],
            "{EHA}_stats/")
            ),
        contigs=os.path.join(
            config["workdir"], 
            "{EHA}_assembly/", 
            "{EHA}_contigs.fasta.gz"
            )
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time='02:00:00'
    message:
        "Creating final assembly summary table for {wildcards.EHA}, uploading files to ERDA"
    shell:
        """
        ### Create the final output summary table
        echo -e "sample\tEHA_number\tEHI_number\tN50\tL50\tnum_contigs\tlargest_contig\tassembly_length\tnum_bins\tassembly_mapping_percent" > {params.stats_dir}/headers.tsv

        #parse QUAST outputs for assembly stats
        for i in {config[workdir]}/bams/*.bam; do 
            cat {params.quast}/{wildcards.EHA}_assembly_report.tsv >> {params.stats_dir}/temp_report.tsv;
        done

        #Add sample IDs
        for ehi in {config[workdir]}/bams/*.bam; do
            echo $(basename ${{ehi/.bam/}}) >> {params.stats_dir}/sample_ids.tsv; 
        done

        for ehi in {config[workdir]}/bams/*.bam; do
            echo {wildcards.EHA} >> {params.stats_dir}/EHA_ids.tsv
        done

        for ehi in {config[workdir]}/bams/*.bam; do
            echo $(basename ${{ehi/_EHA*/}}) >> {params.stats_dir}/EHI_ids.tsv; 
        done
        sed -i 's/PRB.*_//g' {params.stats_dir}/EHI_ids.tsv

        paste {params.stats_dir}/sample_ids.tsv {params.stats_dir}/EHA_ids.tsv {params.stats_dir}/EHI_ids.tsv {params.stats_dir}/temp_report.tsv > {params.stats_dir}/temp2_report.tsv

        for i in {config[workdir]}/bams/*.bam; do 
            cat {params.stats_dir}/{wildcards.EHA}_bins.tsv >> {params.stats_dir}/bins.tsv;
        done

        paste {params.stats_dir}/temp2_report.tsv {params.stats_dir}/bins.tsv > {params.stats_dir}/temp3_report.tsv

        #Grab coverm mapping rate. 'cut -f2' pulls the second column, 'sed -n 3p' prints only the third line (% mapping)
        for i in {input.coverm};
            do cut -f2 $i | sed -n 3p >> {params.stats_dir}/relabun.tsv
        done 

        paste {params.stats_dir}/temp3_report.tsv {params.stats_dir}/relabun.tsv > {params.stats_dir}/temp4_report.tsv

        #Combine them into the final assembly report
        cat {params.stats_dir}/headers.tsv {params.stats_dir}/temp4_report.tsv > {output.stats}


        ### Upload stats to AirTable:
        python {config[codedir]}/airtable/add_asb_stats_airtable.py --report={output.stats} --asb={config[abb]}
        sleep 5

        ### Upload contigs, coverm, & gtdb output to ERDA
        pigz -k -p {threads} {input.contigs}
        
        lftp sftp://erda -e "put {params.contigs} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
        sleep 5
        lftp sftp://erda -e "put {output.stats} -o /EarthHologenomeInitiative/Data/REP/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.gtdb} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"        
        sleep 5
        lftp sftp://erda -e "put {input.tarball} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"

        """