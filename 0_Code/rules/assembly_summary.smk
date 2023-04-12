###############################################################################
## Generate output summary for individual assemblies, update airtable
rule assembly_summary:
    input:
        mw_stats=os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_refinement/",
            "{EHA}_metawrap_70_10_bins.stats",
        ),
        coverm=os.path.join(
            config["workdir"], 
            "{PRB}/",
            "{EHI}/",
            "{EHA}_assembly_coverM.txt"
        ),
    output:
        os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_final_stats.tsv",
        )  
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    params:
        quast=directory(os.path.join(
            config["workdir"], 
            "{PRB}/", 
            "{EHI}/", 
            "{EHA}_QUAST")
        )
    threads: 1
    resources:
        mem_gb=16,
        time='00:05:00'
    message:
        "Creating final coassembly summary table for {wildcards.EHA}"
    shell:
        """
        ### Create the final output summary table
        echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\tassembly_length\tnum_bins\tassembly_mapping_percent" > headers.tsv

        #parse QUAST outputs for assembly stats
        cat {params.quast}/{wildcards.EHA}_assembly_report.tsv >> {wildcards.EHA}_temp_report.tsv

        #Add in the % mapping to assembly stats
        echo {wildcards.EHI}_{wildcards.EHA} >> {wildcards.EHA}_sample_ids.tsv

        paste {wildcards.EHA}_sample_ids.tsv {wildcards.EHA}_temp_report.tsv > {wildcards.EHA}_temp2_report.tsv

        #Add in the # of bins
        cat {wildcards.EHA}_bins.tsv >> {wildcards.EHA}_number_bins.tsv

        paste {wildcards.EHA}_temp2_report.tsv {wildcards.EHA}_number_bins.tsv > {wildcards.EHA}_temp3_report.tsv

        #Grab coverm mapping rate. 'cut -f2' pulls the second column, 'sed -n 3p' prints only the third line (% mapping)
        cut -f2 {input.coverm} | sed -n 3p >> {wildcards.EHA}_relabun.tsv

        paste {wildcards.EHA}_temp3_report.tsv {wildcards.EHA}_relabun.tsv > {wildcards.EHA}_temp4_report.tsv

        #Combine them into the final assembly report
        cat headers.tsv {wildcards.EHA}_temp4_report.tsv > {output}


        ### Upload stats to AirTable:
        python {{config['codedir']}}/airtable/add_asb_stats_airtable.py --report={output} --code={{config['abb']}}
        """