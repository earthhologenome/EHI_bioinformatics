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
        ),
        gtdbtk=os.path.join(
            config["workdir"], 
            config["dmb"] + "_gtdbtk_combined_summary.tsv"
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
    params:
        arch_tree=os.path.join(
            config["workdir"],
            "gtdbtk/classify/gtdbtk.ar53.classify.tree"
        ),
        ct_temp=os.path.join(
            config["workdir"], 
            "count_table_temp.tsv"
        ),
        rm_tax=os.path.join(
            config["workdir"], 
            "rm_tax.tsv"
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
        # Clean up MAGs that don't have GTDBtk classifications
        grep "No bacterial or archaeal marker" {input.gtdbtk} | cut -f1 > {params.rm_tax}
        touch {config["workdir"]}/grep1_done
        grep "Insufficient number of amino acids in MSA" {input.gtdbtk} | cut -f1 >> {params.rm_tax}
        touch {config["workdir"]}/grep2_done
        sed -i 's/.fa//g' {params.rm_tax}
        touch {config["workdir"]}/rm_tax_done
        grep -v -f {params.rm_tax} {input.count_table} > {params.ct_temp}
        touch {config["workdir"]}/ct_gemp_done

        #IF statement, as sometimes we won't have an archaeal tree
        if [ -f {params.arch_tree} ]
        then
            Rscript {config[codedir]}/scripts/prune_tree.R -b {input.tree} -a {params.arch_tree} -i {params.ct_temp} -o {output.tree}

        else
            Rscript {config[codedir]}/scripts/prune_tree_bact.R -b {input.tree} -i {params.ct_temp} -o {output.tree}
        fi

        echo "tree_pruned"

        ## Create list of dereplicated MAG names for updating AirTable
        cut -f1 {params.ct_temp} | sed '1d;' | sed 's/$/.fa/g' > {config[workdir]}/dereplicated_mags.tsv

        ## Remove trailing '.fa'
        sed -i'' 's/\.fa//g' {output.tree}

        ## Remove '.fa' suffix for mag_info, and keep only EHI number for count table
        sed -i'' 's/PRB.....//g' {params.ct_temp}
        sed -i'' 's/_DMB....//g' {params.ct_temp}

        ## Script to split evenness of coverage and counts
        Rscript {config[codedir]}/scripts/split_table.R -i {params.ct_temp} -c {output.counts} -v {output.coverage}

        """