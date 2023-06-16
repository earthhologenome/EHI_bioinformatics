################################################################################
### Run GTDB-tk on refined bins
rule prune_tree:
    input:
        os.path.join(
            config["workdir"],
            config["dmb"] + "_gtdbtk.bac120.classify.tree"
        )
    output:
        os.path.join(
            config["workdir"],
            config["dmb"] + "_pruned.tree"
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
        "Pruning tree"
    shell:
        """
        Rscript {config[codedir]}/scripts/prune_tree.R -i {input} -o {output}
        """