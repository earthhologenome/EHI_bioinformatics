###############################################################################
## Generate output summary for individual assemblies, update airtable
rule assembly_summary:
    input:
        mw_stats=os.path.join(
            config["workdir"],
            "{PRB}_{EHI}_{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        coverm=os.path.join(
            config["workdir"], 
            "coverm/", 
            "{PRB}_{EHI}_{EHA}_assembly_coverM.txt"
            ),
        tarball=os.path.join(
            config["workdir"], 
            "coverm/", 
            "{PRB}_{EHI}_{EHA}_coverm.tar.gz"
            ),
        contigs=os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_assembly/", 
            "{EHA}_contigs.fasta"
            ),
        gtdb=os.path.join(
            config["workdir"], 
            "{PRB}/", 
            "{EHI}/", 
            "{EHA}_gtdbtk_combined_summary.tsv"
            )
    output:
        stats=os.path.join(
            "/projects/ehi/data/REP/",
            "{PRB}_{EHI}_{EHA}_final_stats.tsv",
        ),
        contigs=os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_assembly/", 
            "{EHA}_contigs.fasta.gz"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    params:
        quast=directory(os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_{EHA}_QUAST")
            ),
        stats_dir=directory(os.path.join(
            config["workdir"],
            "{EHA}_stats/")
            )
    threads: 1
    resources:
        load=8,
        mem_gb=16,
        time=estimate_time_summary
    message:
        "Creating final assembly summary table for {wildcards.EHA}, uploading files to ERDA"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt 40 ]
        then
            touch {output.stats}
            touch {output.contigs}

        else

            ### Create the final output summary table
            echo -e "sample\tEHA_number\tEHI_number\tN50\tL50\tnum_contigs\tlargest_contig\tassembly_length\tnum_bins\tassembly_mapping_percent" > {params.stats_dir}/headers.tsv

            #parse QUAST outputs for assembly stats
            cat {params.quast}/{wildcards.EHA}_assembly_report.tsv > {params.stats_dir}/{wildcards.EHA}_temp_report.tsv

            #Add sample IDs
            echo {wildcards.EHI}_{wildcards.EHA} > {params.stats_dir}/{wildcards.EHA}_sample_ids.tsv
            echo {wildcards.EHA} > {params.stats_dir}/{wildcards.EHA}_EHA_ids.tsv
            echo {wildcards.EHI} > {params.stats_dir}/{wildcards.EHA}_EHI_ids.tsv

            paste {params.stats_dir}/{wildcards.EHA}_sample_ids.tsv {params.stats_dir}/{wildcards.EHA}_EHA_ids.tsv {params.stats_dir}/{wildcards.EHA}_EHI_ids.tsv {params.stats_dir}/{wildcards.EHA}_temp_report.tsv > {params.stats_dir}/{wildcards.EHA}_temp2_report.tsv

            paste {params.stats_dir}/{wildcards.EHA}_temp2_report.tsv {params.stats_dir}/{wildcards.EHA}_bins.tsv > {params.stats_dir}/{wildcards.EHA}_temp3_report.tsv

            #Grab coverm mapping rate. 'cut -f2' pulls the second column, 'sed -n 3p' prints only the third line (% mapping)
            cut -f2 {input.coverm} | sed -n 3p > {params.stats_dir}/{wildcards.EHA}_relabun.tsv

            paste {params.stats_dir}/{wildcards.EHA}_temp3_report.tsv {params.stats_dir}/{wildcards.EHA}_relabun.tsv > {params.stats_dir}/{wildcards.EHA}_temp4_report.tsv

            #Combine them into the final assembly report
            cat {params.stats_dir}/headers.tsv {params.stats_dir}/{wildcards.EHA}_temp4_report.tsv > {output.stats}


            ### Upload stats to AirTable:
            python {config[codedir]}/airtable/add_asb_stats_airtable.py --report={output.stats} --asb={config[abb]}
            sleep 5

            ### Upload contigs, coverm, & gtdb output to ERDA
            pigz -k -p {threads} {input.contigs}
            
            lftp sftp://erda -e "put {output.contigs} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
            sleep 5
            lftp sftp://erda -e "put {output.stats} -o /EarthHologenomeInitiative/Data/REP/; bye"
            sleep 5
            lftp sftp://erda -e "put {input.gtdb} -o /EarthHologenomeInitiative/Data/REP/; bye"
            sleep 5
            lftp sftp://erda -e "put {input.tarball} -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
            
        fi
        """