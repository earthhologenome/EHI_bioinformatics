################################################################################
### Run GTDB-tk on refined bins
rule gtdbtk:
    input:
        stats=os.path.join(
            config["workdir"],
            "{PRB}",
            "{EHI}",
            "{EHA}_refinement",
            "{EHA}_metawrap_70_10_bins.stats",
        )
    output:
        os.path.join(config["workdir"], "{PRB}", "{EHI}", "{EHA}", "gtdbtk/classify/gtdbtk.bac120.summary.tsv")
    params:
        GTDB_data=expand("{GTDB_data}", GTDB_data=config['GTDB_data']),
        outdir="{config['workdir'']}/{PRB}/{EHI}/{EHA}/gtdbtk",
        bins="{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/metawrap_70_10_bins"
    conda:
        f"{config['codedir']}/conda_envs/GTDB-tk.yaml"
    threads:
        16
    resources:
        mem_gb=96,
        time='02:00:00'
    benchmark:
        "{{config['logdir']}}/gtdb-tk_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/gtdb-tk_log_{PRB}_{EHI}_{EHA}.log"
    message:
        "Annotating taxonomy to {wildcards.EHA}'s bins using GTDB-tk"
    shell:
        """
        # Specify path to reference data:
        export GTDBTK_DATA_PATH={params.GTDB_data}

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.bins} \
        --extension "gz" \
        --out_dir {params.outdir} \
        --cpus {threads}

        # Create a merged summary output for DRAM:
        if [ -s "{params.outdir}/classify/gtdbtk.ar122.summary.tsv" ]
        then
        sed '1d;' {params.outdir}/classify/gtdbtk.ar122.summary.tsv > {params.outdir}/ar122.tsv
        cat {output} {params.outdir}/ar122.tsv > {params.outdir}/gtdbtk_combined_summary.tsv
        rm {params.outdir}/ar122.tsv

        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {output} > {params.outdir}/gtdbtk_combined_summary.tsv
        fi
        """