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
        ),
        taxonomy=os.path.join(
            config["workdir"],
            config["dmb"] + "_taxon_table.tsv"
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
        mem_gb=348,
        time='08:00:00'
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

        # Split taxonomy into taxon group columns:
        cut -f1 {output.combined} > {params.outdir}/mag_names.tsv
        cut -f2 {output.combined} | sed '1d;' | tr ';' '\t' > {params.outdir}/taxonomy.tsv
        echo -e 'mag_name\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies' > {params.outdir}/gtdb_headers.tsv
        paste {params.outdir}/mag_names.tsv {params.outdir}/taxonomy.tsv > {params.outdir}/gtdb_temp.tsv
        cat {params.outdir}/gtdb_headers.tsv {params.outdir}/gtdb_temp.tsv > {output.taxonomy}

        cp {params.outdir}/classify/gtdbtk.bac120.classify.tree {output.tree}
        """