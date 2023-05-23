################################################################################
### Run GTDB-tk on refined bins
rule gtdbtk:
    input:
        stats=os.path.join(
            config["workdir"],
            "{PRB}_{EHI}_{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
            )
    output:
        bac=os.path.join(
            config["workdir"], 
            "{PRB}/", 
            "{EHI}/", 
            "{EHA}/", 
            "gtdbtk/classify/gtdbtk.bac120.summary.tsv"
            ),
        combined=os.path.join(
            config["workdir"], 
            "{PRB}/", 
            "{EHI}/",
            "{EHA}_gtdbtk_combined_summary.tsv"
            )
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
        mem_gb=72,
        time=estimate_time_gtdb
    benchmark:
        os.path.join(config["logdir"] + "/gtdb-tk_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/gtdb-tk_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Annotating taxonomy to {wildcards.EHA}'s bins using GTDB-tk"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt 40 ]
        then
            touch {output.bac}
            touch {output.combined}

        else

            # Create temp folder
            export TMPDIR={config[workdir]}/tmpdir
            mkdir -p $TMPDIR
            
            # Specify path to reference data:
            export GTDBTK_DATA_PATH={params.GTDB_data}

            # Get version number for AirTable:
            gtdbtk --version | cut -f3 -d ' ' > {params.outdir}/version.tsv

            for i in {params.bins}/*.fa.gz;
                do cat {params.outdir}/version.tsv >> {params.outdir}/version_catted.tsv;
            done

            # Run GTDB-tk:
            gtdbtk classify_wf \
            --genome_dir {params.bins} \
            --extension "gz" \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --skip_ani_screen

            # Create a merged summary output for DRAM:
            if [ -s "{params.outdir}/classify/gtdbtk.ar53.summary.tsv" ]
            then
                sed '1d;' {params.outdir}/classify/gtdbtk.ar53.summary.tsv > {params.outdir}/ar53.tsv
                cat {output.bac} {params.outdir}/ar53.tsv > {output.combined}
            rm {params.outdir}/ar53.tsv

            # Otherwise, just use the bacterial summary (if no archaeal bins)
            else
                cat {output.bac} > {output.combined}
            fi 
            
            # Parse the gtdb output for uploading to the EHI MAG database
            cut -f2 {output.combined} | sed '1d;' | tr ';' '\t' > {params.outdir}/taxonomy.tsv
            cut -f1,6,11,12 {output.combined} | sed '1d;' > {params.outdir}/id_ani.tsv
            sed -i 's@N/A@0@g' {params.outdir}/id_ani.tsv
            echo -e 'mag_name\tfastani_ani\tclosest_placement_ani\tclosest_placement_af\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > {params.outdir}/gtdb_headers.tsv
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

            paste {params.outdir}/mag_names.tsv {params.outdir}/n_contigs.tsv {params.outdir}/EHA.tsv {params.outdir}/version_catted.tsv > {params.outdir}/ncontigs_temp.tsv
            echo -e 'mag_id\tcontigs\teha_number\tGTDB_version' > {params.outdir}/ncontigs_header.tsv
            cat {params.outdir}/ncontigs_header.tsv {params.outdir}/ncontigs_temp.tsv > {params.outdir}/ncontigs_temp2.tsv

            cut -f1,2,3,4,6,7 {params.refinement}/{wildcards.EHA}_metawrap_50_10_bins.stats > {params.outdir}/mw_stats.tsv

            #combine into final table for upload to airtable:
            paste {params.outdir}/gtdb_airtable.tsv {params.outdir}/ncontigs_temp2.tsv {params.outdir}/mw_stats.tsv > /projects/ehi/data/REP/{config[abb]}_{wildcards.EHA}_mags.tsv

            #update EHI MAG airtable with stats:
            python {config[codedir]}/airtable/add_mag_stats_airtable.py --report=/projects/ehi/data/REP/{config[abb]}_{wildcards.EHA}_mags.tsv

        fi
        """