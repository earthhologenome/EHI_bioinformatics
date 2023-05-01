################################################################################
### Run GTDB-tk on refined bins
rule gtdbtk:
    input:
        stats=os.path.join(
            config["workdir"],
            "{PRB}_{EHI}_{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            )
    output:
        os.path.join(
            config["workdir"], 
            "{PRB}/", 
            "{EHI}/", 
            "{EHA}/", 
            "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    params:
        GTDB_data=expand("{GTDB_data}", GTDB_data=config['GTDB_data']),
        outdir=os.path.join(config["workdir"] + "/{PRB}" + "/{EHI}" + "/{EHA}" + "/gtdbtk"),
        refinement=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_refinement"),
        bins=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_refinement" + "/metawrap_50_10_bins")
    conda:
        f"{config['codedir']}/conda_envs/GTDB-tk.yaml"
    threads:
        16
    resources:
        mem_gb=96,
        time='08:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/gtdb-tk_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/gtdb-tk_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Annotating taxonomy to {wildcards.EHA}'s bins using GTDB-tk"
    shell:
        """
        # Specify path to reference data:
        export GTDBTK_DATA_PATH={params.GTDB_data}

        # Get version number for AirTable:
        gtdbtk --version | cut -f3 -d ' ' > {params.outdir}/version.tsv

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.bins} \
        --extension "gz" \
        --out_dir {params.outdir} \
        --cpus {threads} \
        --mash_db /projects/ehi/data/0_Environments/databases/{wildcards.EHA}_gtdb-tk-r207.msh

        # Create a merged summary output for DRAM:
        if [ -s "{params.outdir}/classify/gtdbtk.ar122.summary.tsv" ]
        then
        sed '1d;' {params.outdir}/classify/gtdbtk.ar122.summary.tsv > {params.outdir}/ar122.tsv
        cat {output} {params.outdir}/ar122.tsv > {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv
        rm {params.outdir}/ar122.tsv

        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {output} > {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv
        fi

        # Parse the gtdb output for uploading to the EHI MAG database
        cut -f2 {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv | sed '1d;' | tr ';' '\t' > {params.outdir}/taxonomy.tsv
        cut -f1,11 {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv | sed '1d;' > {params.outdir}/id_ani.tsv
        echo -e 'mag_name\tclosest_placement_ani\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > {params.outdir}/gtdb_headers.tsv
        paste {params.outdir}/id_ani.tsv {params.outdir}/taxonomy.tsv > {params.outdir}/gtdb_temp.tsv
        cat {params.outdir}/gtdb_headers.tsv {params.outdir}/gtdb_temp.tsv > {params.outdir}/gtdb_airtable.tsv

        # Get the # contigs per MAG, also completeness/contamination/size from metawrap stats
        for i in {params.bins}/*.fa.gz;
            do echo $(basename $i) >> {params.outdir}/mag_names.tsv && zcat $i | grep '>' | wc -l >> {params.outdir}/n_contigs.tsv;
        done

        # Get the EHA number
        for i in {params.bins}/*.fa.gz;
            do echo {wildcards.EHA} >> {params.outdir}/EHA.tsv;
        done

        paste {params.outdir}/mag_names.tsv {params.outdir}/n_contigs.tsv {params.outdir}/EHA.tsv {params.outdir}/version.tsv > {params.outdir}/ncontigs_temp.tsv
        echo -e 'mag_id\tcontigs\teha_number\tGTDB_version' > {params.outdir}/ncontigs_header.tsv
        cat {params.outdir}/ncontigs_header.tsv {params.outdir}/ncontigs_temp.tsv > {params.outdir}/ncontigs_temp2.tsv

        cut -f1,2,3,4,6,7 {params.refinement}/{wildcards.EHA}_metawrap_50_10_bins.stats > {params.outdir}/mw_stats.tsv

        #combine into final table for upload to airtable:
        paste {params.outdir}/gtdb_airtable.tsv {params.outdir}/ncontigs_temp2.tsv {params.outdir}/mw_stats.tsv > {params.outdir}/airtable_mag.tsv

        #update EHI MAG airtable with stats:
        python {config[codedir]}/airtable/add_mag_stats_airtable.py --report={params.outdir}/airtable_mag.tsv
        """