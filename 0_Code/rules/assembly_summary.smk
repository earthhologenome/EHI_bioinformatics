###############################################################################
## Generate output summary for individual assemblies, update airtable
rule assembly_summary:
    input:
        mw_stats=os.path.join(
            config["workdir"],
            "{PRB}",
            "{EHI}",
            "{EHA}_refinement",
            "{EHA}_metawrap_70_10_bins.stats",
        ),
        coverm=os.path.join(
            config["workdir"], 
            "{PRB}",
            "{EHI}",
            "{EHA}_assembly_coverM.txt"
        ),
    output:
        os.path.join(
            config["workdir"],
            "{PRB}",
            "{EHI}",
            "{EHA}_final_stats.tsv",
        )  
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    params:
        quast=directory(os.path.join(
            config["workdir"], 
            "{PRB}", 
            "{EHI}", 
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
        #Create the final output summary table
        #parse QUAST outputs for assembly stats
        echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\tassembly_mapping_percent" > headers.tsv


        for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
            do cat 3_Outputs/2_Coassemblies/{params.group}_QUAST/{params.group}_assembly_report.tsv >> {params.group}_temp_report.tsv;
        done


        #Add in the % mapping to assembly stats
        for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
            do echo $(basename ${{sample/.bam/}}) >> {params.group}_sample_ids.tsv;
        done

        paste {params.group}_sample_ids.tsv {params.group}_temp_report.tsv > {params.group}_temp2_report.tsv

        #Add in the # of bins
        for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
            do cat {params.group}_bins.tsv >> {params.group}_number_bins.tsv;
        done

        paste {params.group}_temp2_report.tsv {params.group}_number_bins.tsv > {params.group}_temp3_report.tsv

        ls -l 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam | wc -l > {params.group}_n_samples.tsv

        nsamples=$( cat {params.group}_n_samples.tsv )
        nsamples1=$(( nsamples + 1 ))
        echo $nsamples1
        for sample in `seq 2 $nsamples1`;
            do cut -f"$sample" 3_Outputs/6_Coassembly_CoverM/{params.group}_assembly_coverM.txt | sed -n 3p >> {params.group}_relabun.tsv;
        done

        paste {params.group}_temp3_report.tsv {params.group}_relabun.tsv > {params.group}_temp4_report.tsv
        #Combine them into the final assembly report
        cat headers.tsv {params.group}_temp4_report.tsv > {output}


        # Upload stats to AirTable:


        # Log on the AirTable that the pipeline has finished:


        """
