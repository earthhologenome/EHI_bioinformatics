################################################################################
### Run GTDB-tk on refined bins
rule gtdbtk_full_tree:
    input:
        os.path.join(
            config["workdir"],
            "drep/",
            "figures/",
            config["dmb"] + "_Primary_clustering_dendrogram.pdf"
        )
    output:
        bac=os.path.join(
            config["workdir"], 
            "gtdbtk/classify/gtdbtk.bac120.summary.tsv"
        ),
        combined=os.path.join(
            config["workdir"], 
            config["dmb"] + "_gtdbtk_combined_summary.tsv"
        ),
        tree=os.path.join(
            config["workdir"],
            config["dmb"] + "_gtdbtk.bac120.classify.tree"
        )
    params:
        GTDB_data=expand("{GTDB_data}", GTDB_data=config['GTDB_data']),
        outdir=os.path.join(config["workdir"] + "/gtdbtk"),
        dereplicated_mags=os.path.join(config["workdir"] + "/drep/dereplicated_genomes"),
    conda:
        f"{config['codedir']}/conda_envs/GTDB-tk.yaml"
    threads:
        16
    resources:
        mem_gb=456,
        time=estimate_time_gtdb
    benchmark:
        os.path.join(config["logdir"] + "/gtdb-tk_benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/gtdb-tk_log.log")
    message:
        "Creating phylogenetic tree from derepliated mags using GTDB-tk"
    shell:
        """
        # Specify path to reference data:
        export GTDBTK_DATA_PATH={params.GTDB_data}

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.dereplicated_mags} \
        --extension "gz" \
        --out_dir {params.outdir} \
        --cpus {threads} \
        --skip_ani_screen \
        --full_tree

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

        cp {params.outdir}/classify/gtdbtk.bac120.classify.tree {output.tree}
        """