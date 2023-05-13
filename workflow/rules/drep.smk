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
        # Load dRep module, as the conda recipe is cooked atm
        module load drep/3.4.0

        # Massage genome info file:
        sed 's/.fa/.fa.gz/g' mags.csv > mags_formatted.csv

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
        """