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
        # remove folder (in case of restart/server interruption)
        rm -rf {config[workdir]}/drep

        # mjolnir drep module needs to have dependencies installed, here's a temp solution
        module load conda/25.1.1
        module load fastani/1.33
        source activate /projects/ehi/data/0_Environments/conda/drep

        # Massage genome info file:
        cut -f1,3,4 -d ',' mags.csv | sed 's/.fa/.fa.gz/g' > mags_formatted.csv

        # Dereplicate these suckers:
        dRep dereplicate \
                {config[workdir]}/drep \
                -p {threads} \
                -comp 50 \
                -sa {config[ani]} \
                -g {config[magdir]}/*.fa.gz \
                --genomeInfo mags_formatted.csv
                2> {log}

        for i in {config[workdir]}/drep/figures/*;
            do mv $i {config[workdir]}/drep/figures/{config[dmb]}_$(basename "$i");
        done

        for i in {config[workdir]}/drep/data_tables/*;
            do mv $i {config[workdir]}/drep/data_tables/{config[dmb]}_$(basename "$i");
        done

        tar -cvzf {config[workdir]}/{config[dmb]}_drep.tar.gz {config[workdir]}/drep/data_tables {config[workdir]}/drep/figures

        """