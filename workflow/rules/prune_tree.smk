################################################################################
### Run GTDB-tk on refined bins
rule prune_tree:
    input:
        tree=os.path.join(
            config["workdir"],
            config["dmb"] + "_gtdbtk.bac120.classify.tree"
        ),
        count_table=os.path.join(
            config["workdir"], 
            "coverm/", 
            config["dmb"] + "_count_table.tsv"
        )
    output:
        tree=os.path.join(
            config["workdir"],
            config["dmb"] + ".tree"
        ),
        counts=os.path.join(
            config["workdir"],
            config["dmb"] + "_counts.tsv"
        ),
        coverage=os.path.join(
            config["workdir"],
            config["dmb"] + "_coverage.tsv"
        )
    conda:
        f"{config['codedir']}/conda_envs/R_tidyverse.yaml"
    threads:
        2
    resources:
        mem_gb=32,
        time='00:10:00'
    log:
        os.path.join(config["logdir"] + "/gtdb-tk_log.log")
    message:
        "Pruning tree & splitting count table"
    shell:
        """
        Rscript {config[codedir]}/scripts/prune_tree.R -i {input.tree} -o {output.tree}

        ## Remove trailing '.fa'
        sed -i'' 's/\.fa//g' {output.tree}

        ## Remove '.fa' suffix for mag_info, and keep only EHI number for count table
        sed -i'' 's/PRB.....//g' {input.count_table}
        sed -i'' 's/_DMB....//g' {input.count_table}

        ## Script to split evenness of coverage and counts
        Rscript {config[codedir]}/scripts/split_table.R -i {input.count_table} -c {output.counts} -v {output.coverage}

        """