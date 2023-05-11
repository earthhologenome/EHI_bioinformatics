################################################################################
### Dereplicate MAGs using dRep
rule drep:
    input:
        downloaded=os.path.join(
            config["magdir"],
            "mags_downloaded"
        )
    output:
        os.path.join(
            config["workdir"],
            "drep/",
            "figures/",
            config["dmb"] + "_Primary_clustering_dendrogram.pdf"
        )
    conda:
        f"{config['codedir']}/conda_envs/drep.yaml"
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='08:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/drep.benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/drep.log")
    message:
        "Dereplicating MAGs with dRep"
    shell:
        """
        # Massage genome info file:
        sed -i 's/.fa/.fa.gz/g' mags.csv

        # Dereplicate these suckers:
        dRep dereplicate \
                {config[workdir]}/drep \
                -p {threads} \
                -comp 50 \
                -sa {config[ani]} \
                -g {config[magdir]}/*.fa.gz \
                --genomeInfo mags.csv
                2> {log}

        for i in {config[workdir]}/drep/figures/*;
            do mv $i {config[dmb]}_"$i";
        done
        """